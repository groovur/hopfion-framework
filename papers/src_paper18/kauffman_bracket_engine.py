#!/usr/bin/env python3
"""
kauffman_bracket_engine.py
===================================================================
Minimal, from-first-principles Kauffman-bracket / Temperley-Lieb (TL_n)
calculator, built to independently verify Jones-polynomial claims about
the torus knot/link tower (T(2,n)) used throughout Papers XVIII/XIX,
without relying on memorized R-matrix / quantum-group normalization
conventions (which are easy to get subtly wrong by hand -- see
jones_T23_satellite_check.py for a worked example of exactly that).

Method
------
A tangle diagram on n strands (top points 0..n-1, bottom points
n..2n-1) is represented as a perfect non-crossing matching of the 2n
boundary points. A general tangle is a formal linear combination of
such matchings with coefficients in Z[A, A^{-1}] (sympy).

  - identity_tangle(n)      : n parallel vertical strands
  - e_tangle(n, i)           : Temperley-Lieb generator e_i (cap/cup
                               connecting strands i, i+1)
  - crossing_tangle(n, i, +-1): Kauffman skein relation for a single
                               crossing between adjacent strands i,i+1:
                                 positive crossing = A*identity + A^-1*e_i
                                 negative crossing = A^-1*identity + A*e_i
  - compose_tangles(t1, t2, n): stack t1 on top of t2 (vertical
                               composition), correctly counting any
                               closed loops formed at the gluing
                               interface (each contributes a factor
                               delta = -A^2 - A^-2)
  - close_up_bracket(tangle, n): trace closure (connect top i to
                               bottom i for all i), returns the
                               Kauffman bracket normalized so that a
                               single unknot = 1.

Jones polynomial from the bracket: for a diagram with writhe w,
    V(t) = [ (-A^3)^{-w} * <bracket> ]  with A = t^{-1/4}.

Validated against the two smallest nontrivial torus links (see
jones_T23_satellite_check.py): the trefoil T(2,3) and the Hopf link
T(2,2), both reproduced exactly against their standard tabulated
Jones polynomials before being trusted on the new (uncomputed-by-the-
papers) satellite diagram.
"""

import sympy as sp

A = sp.symbols('A')
delta = -A**2 - A**-2

# ---------------------------------------------------------------
# Temperley-Lieb tangles as linear combinations of non-crossing
# matchings of 2n boundary points: top = 0..n-1, bottom = n..2n-1.
# A "tangle" is a dict: {matching(frozenset of frozenset pairs): coeff}
# ---------------------------------------------------------------

def identity_diagram(n):
    return frozenset(frozenset((i, i + n)) for i in range(n))

def e_gen_diagram(n, i):
    pairs = [frozenset((i, i + 1)), frozenset((i + n, i + 1 + n))]
    for j in range(n):
        if j != i and j != i + 1:
            pairs.append(frozenset((j, j + n)))
    return frozenset(pairs)

def add_tangles(*ts):
    out = {}
    for t in ts:
        for k, v in t.items():
            out[k] = out.get(k, 0) + v
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}

def scale_tangle(t, s):
    return {k: sp.expand(v * s) for k, v in t.items()}

def compose_diagrams(d1, d2, n):
    """Stack d1 on top of d2. d1: top 0..n-1, bottom n..2n-1 (glue labels).
       d2: top 0..n-1 (glue labels), bottom n..2n-1.
       Returns (result_matching, num_closed_loops)."""
    def tag1(x):
        return ('T', x) if x < n else ('G', x - n)
    def tag2(x):
        return ('G', x) if x < n else ('B', x - n)

    adj = {}
    def add_edge(u, v):
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)

    for pair in d1:
        a, b = tuple(pair)
        add_edge(tag1(a), tag1(b))
    for pair in d2:
        a, b = tuple(pair)
        add_edge(tag2(a), tag2(b))

    visited_nodes = set()
    loops = 0
    result_pairs = []

    all_nodes = list(adj.keys())
    for start in all_nodes:
        if start in visited_nodes:
            continue
        kind, _ = start
        if kind == 'G':
            continue  # only start traces from boundary nodes; G-only cycles handled after
        visited_nodes.add(start)
        prev = None
        cur = start
        while True:
            nbrs = adj[cur]
            cands = [x for x in nbrs if x != prev] if prev is not None else list(nbrs)
            nxt = cands[0] if cands else prev
            if nxt[0] != 'G':
                result_pairs.append(frozenset((start, nxt)))
                visited_nodes.add(nxt)
                break
            visited_nodes.add(nxt)
            prev, cur = cur, nxt

    for start in all_nodes:
        if start in visited_nodes or start[0] != 'G':
            continue
        cur = start
        prev = None
        while True:
            visited_nodes.add(cur)
            nbrs = adj[cur]
            cands = [x for x in nbrs if x != prev] if prev is not None else list(nbrs)
            nxt = cands[0] if cands else nbrs[0]
            if nxt == start:
                break
            prev, cur = cur, nxt
        loops += 1

    final_pairs = []
    for pair in result_pairs:
        u, v = tuple(pair)
        def relabel(node):
            k, x = node
            return x if k == 'T' else x + n
        final_pairs.append(frozenset((relabel(u), relabel(v))))

    return frozenset(final_pairs), loops

def compose_tangles(t1, t2, n):
    out = {}
    for d1, c1 in t1.items():
        for d2, c2 in t2.items():
            d, loops = compose_diagrams(d1, d2, n)
            coeff = sp.expand(c1 * c2 * delta**loops)
            out[d] = out.get(d, 0) + coeff
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}

def identity_tangle(n):
    return {identity_diagram(n): sp.Integer(1)}

def e_tangle(n, i):
    return {e_gen_diagram(n, i): sp.Integer(1)}

def crossing_tangle(n, i, sign=+1):
    idt = identity_tangle(n)
    e = e_tangle(n, i)
    if sign == +1:
        return add_tangles(scale_tangle(idt, A), scale_tangle(e, A**-1))
    else:
        return add_tangles(scale_tangle(idt, A**-1), scale_tangle(e, A))

def close_up_bracket(tangle, n):
    """Trace closure: connect top i to bottom (n+i) for all i, count loops, return bracket
       normalized so a single unknot = 1 (i.e. delta^(loops-1))."""
    total = 0
    for d, c in tangle.items():
        adj = {}
        def add_edge(u, v):
            adj.setdefault(u, []).append(v)
            adj.setdefault(v, []).append(u)
        for pair in d:
            a, b = tuple(pair)
            add_edge(a, b)
        for i in range(n):
            add_edge(i, n + i)
        visited = set()
        loops = 0
        for node in list(adj.keys()):
            if node in visited:
                continue
            cur = node
            prev = None
            while True:
                visited.add(cur)
                nbrs = adj[cur]
                cands = [x for x in nbrs if x != prev] if prev is not None else list(nbrs)
                nxt = cands[0] if cands else nbrs[0]
                if nxt == node:
                    break
                prev, cur = cur, nxt
            loops += 1
        total += c * delta**(loops - 1)
    return sp.simplify(total)

if __name__ == "__main__":
    # sanity check: TL_2 relation e*e = delta*e
    n = 2
    e = e_tangle(n, 0)
    ee = compose_tangles(e, e, n)
    print("e*e should be delta*e:", ee)
