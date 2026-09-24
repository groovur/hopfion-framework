#!/usr/bin/env python3
r"""
Do the E_8 condensate-scale mass ratios (pred/PDG) correlate with the twist chirality,
or with something else? (2026-09-03). See notes/quark_generation_e8_ribbon_twist.md.

ratio = E_8 prediction / PDG at the condensate scale (from e8_strange_anchor.py).
Tw = m-15 (twist). isospin = up/down.
"""
import statistics as st

# q : (Tw, isospin, ratio=pred/PDG, PDG mass MeV)
D = {
 'u': (+4,'up',   0.674,   2.16),
 'c': (-4,'up',   0.335,   1275),
 't': (-14,'up',  1.935,   172570),
 'd': (+2,'down', 1.260,   4.67),
 's': (-2,'down', 1.131,   93.4),
 'b': (-8,'down', 1.211,   4180),
}

print("="*70); print("E_8 mass ratios vs twist chirality vs isospin"); print("="*70)
print("\n q  Tw   isospin   pred/PDG")
for q,(tw,iso,r,m) in sorted(D.items(), key=lambda kv: kv[1][0]):
    print(f" {q:>2} {tw:>+3}   {iso:>4}     {r:5.3f}")

up   = [D[q][2] for q in ('u','c','t')]
down = [D[q][2] for q in ('d','s','b')]
print("\n[1] By ISOSPIN:")
print(f"    down-type {['%.3f'%r for r in down]}  mean={st.mean(down):.3f}  spread(max/min)={max(down)/min(down):.2f}")
print(f"    up-type   {['%.3f'%r for r in up]}  mean={st.mean(up):.3f}  spread(max/min)={max(up)/min(up):.2f}")
print("    => DOWN-type clusters tight (~1.20, +-5%): the down mass RATIOS d:s:b are correct up to a")
print("       single uniform ~1.2 factor (a scale/running normalisation).")
print("       UP-type scatters over a factor ~6: the up ratios u:c:t are WRONG, not a uniform factor.")

def pearson(xs, ys):
    n=len(xs); mx=st.mean(xs); my=st.mean(ys)
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return cov/(sx*sy)

qs = list(D)
tw_all = [D[q][0] for q in qs]; r_all = [D[q][2] for q in qs]
absTw  = [abs(t) for t in tw_all]
iso_num = [1 if D[q][1]=='down' else 0 for q in qs]
print("\n[2] Pearson correlation of ratio with (all 6 points):")
print(f"    signed Tw : {pearson(tw_all, r_all):+.2f}   |Tw| : {pearson(absTw, r_all):+.2f}   isospin(down=1) : {pearson(iso_num, r_all):+.2f}")

print("\n[2b] ROBUSTNESS -- leave-one-out (is the |Tw| correlation real or one-point?):")
for drop in ('t','c'):
    ks=[q for q in qs if q!=drop]
    a=[abs(D[q][0]) for q in ks]; r=[D[q][2] for q in ks]
    print(f"    drop {drop}: |Tw| Pearson = {pearson(a,r):+.2f}")
print("    => the +0.67 |Tw| correlation is driven ENTIRELY by top (drop t -> collapses). |Tw|=4 (c,u)")
print("       dips BELOW |Tw|=2 (s,d), so it is non-monotone: NOT a robust twist law. c is an outlier too.")

print("\n[3] The ROBUST feature is isospin VARIANCE, not a Pearson trend:")
print(f"    down {['%.2f'%D[q][2] for q in ('d','s','b')]} spread x{max(down)/min(down):.2f}  (TIGHT ~1.20)")
print(f"    up   {['%.2f'%D[q][2] for q in ('u','c','t')]} spread x{max(up)/min(up):.2f}  (SCATTERED)")
print("    All 3 down cluster; all 3 up scatter -> robust. Pearson-of-mean misses it because the MEANS")
print("    are close (1.20 vs 0.98); the signal is the VARIANCE. => down ratios correct up to a uniform")
print("    ~1.2 (running); up ratios genuinely wrong.")

print("\n[4] WHY, and twist's real role:")
print("    Anchor is DOWN (strange, n_s=0 = muon) -> down sector calibrated -> internal ratios right.")
print("    UP has no anchor (an up-anchor needs a ~meV neutrino resonance; impossible) -> extrapolated,")
print("    scattered; and the up hierarchy is steeper (m_t/m_u~1e5 vs m_b/m_d~1e3), the harder sector.")
print("    Twist's role is INDIRECT: twist -> isospin (strange down via muon) -> which sector is anchored.")
print("    HONEST: no direct twist-chirality<->ratio law (the apparent Tw correlation is top-driven);")
print("    the real, robust split is isospin, set by the anchor.")
print("="*70)
