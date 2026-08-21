"""Two follow-up structural checks on the p^2/p^4 fingerprint null.

A. Shell-dependence: the number of Russell-Saunders terms of l^2 is 2l+1
   (s:1, p:3, d:5, f:7), triplet count = l. So the p-shell "3 terms = 3
   SU(2)_3 non-vacuum primaries" coincidence is 2l+1 evaluated at l=1; it
   does not recur (d->5, f->7), confirming SO(3) origin, not a fixed
   condensate integer.
B. The WZW<->SO(3) relation: the SU(2)_3 quantum dimension
   [n]_q = sin(n pi/(k+2))/sin(pi/(k+2)) at q=exp(i pi/5) gives
   [2]_q=[3]_q=phi, [4]_q=1, [5]_q=0 (the j<=3/2 truncation), while the
   classical limit q->1 gives [n]_1 = n, the ordinary SO(3) dimension
   2l+1. The condensate's phi is the q-deformation of the SO(3) dim 2;
   atoms sit at q=1 (integers), the condensate at q=exp(i pi/5) (phi +
   truncation). Same algebra, same formula, different q -- which is why
   atomic multiplets show 2l+1 and never phi.
"""
import sympy as sp

def terms_l2(l):
    return [(L, 0 if (2*l-L)%2==0 else 1) for L in range(2*l, -1, -1)]

def qint(n, K=5):
    return sp.simplify(sp.sin(n*sp.pi/K)/sp.sin(sp.pi/K))

if __name__ == '__main__':
    for name,l in [('s2',0),('p2',1),('d2',2),('f2',3)]:
        T=terms_l2(l)
        assert len(T)==2*l+1
        print(f"{name}: {len(T)} terms (=2l+1), {sum(1 for _,S in T if S==1)} triplets")
    phi=(1+sp.sqrt(5))/2
    assert sp.simplify(qint(2)-phi)==0 and sp.simplify(qint(3)-phi)==0
    assert sp.simplify(qint(4)-1)==0 and sp.simplify(qint(5))==0
    print("q-deform verified: [2]=[3]=phi, [4]=1, [5]=0 (truncation); [n]_{q=1}=n (SO(3)).")
