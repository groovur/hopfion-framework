# delta1 numeric search results

Engine: mpmath 1.3.0 (Python 3.11, via `python3.11`), mp.dps=60. numpy 2.2.6 used for
vectorized Test B search and Test A calibration.

## Definitions

D1 = 16*phi**8/5 - 360/phi**2 = 12.8241199939550791529994956410006038901160816130218947...
D1_alt = (2136*phi - 3392)/5 = 12.8241199939550791529994956410006038901160816130218947...
|D1 - D1_alt| = 3.73e-60 -> AGREE to 50+ digits. VERIFIED, proceeded.

C = 4*pi + 15/8 - phi = 12.82333662560927810564598669875237341906836841769466042...

delta = D1 - C (50 digits) =
0.00078336834580104735350894224823047104771319532723428

## Engine validation

Planted T0 = 3*pi/2 + 7/8 - 2*phi -> (a,b,c) = (7/8, 3/2, -2), expected total complexity
= comp(7/8)+comp(3/2)+comp(-2) = 15+5+3 = 23.
Enumerated all expressions up to complexity 23 (399,575 distinct triples): target found
exactly at complexity 23. PASS.

## Test A (grammar E = a + b*pi + c*phi, comp<=30)

Enumerated 1,936,179 distinct (a,b,c) triples with total complexity <= 30 (10.7s enumerate,
3.0s evaluate in double precision).

Hits with relerr <= 6.15e-5 (candidate's own tolerance): 40 total, sorted by complexity
(full list in testA_hits.txt). Candidate (a=15/8, b=4, c=-1), complexity 30, relerr=6.109e-05
IS present (as required).

Complexity-30 hits (11 of the 40):
```
comp=30  a=  -3/2 b=  -1/3 c=  19/2  relerr=4.170e-07
comp=30  a=   1/4 b=  11/5 c=   7/2  relerr=3.921e-05
comp=30  a=   7/2 b=     3 c= -1/16  relerr=3.658e-05
comp=30  a=    10 b=   5/9 c=   2/3  relerr=7.908e-06
comp=30  a=   7/4 b=   8/3 c=   5/3  relerr=1.433e-05
comp=30  a= -11/2 b=     3 c=  11/2  relerr=1.209e-05
comp=30  a=  -9/4 b=     6 c=  -7/3  relerr=1.816e-06
comp=30  a=   6/7 b=   7/2 c=   3/5  relerr=4.542e-05
comp=30  a=  6/11 b=   4/3 c=     5  relerr=2.298e-05
comp=30  a= -3/14 b=    -1 c=    10  relerr=2.663e-05
comp=30  a=  15/8 b=     4 c=    -1  relerr=6.109e-05   <-- candidate C
```
(Full 40-row list, all complexities 20-30, in testA_hits.txt.)

Calibration: 200 pseudo-random targets uniform in [11.5,14.0], seed=42, same tolerance
6.15e-5, counted against the same sorted array of 1,936,179 enumerated values:
- mean hits = 35.74
- min = 20, median = 35.0, max = 56

Verdict data: D1 observed 40 hits vs expected-by-chance mean 35.74 (range 20-56 over 200
draws). 40 is within one distribution-width of the mean for an arbitrary target in this
range -- not a statistically surprising hit count.

## Test B (basis {pi,phi,sqrt5,pi^2,phi^2,pi*phi,pi/phi,1/pi,1/phi,log(phi),log(2)},
q0 + q1*B_i + q2*B_j, den<=12, |num|<=48, total comp<=40)

412,179 candidate expressions evaluated (1,131 single-basis-term + 411,048 two-basis-term,
after budget pruning); double-precision search, frontier re-verified at 60-digit mpmath.

Pareto frontier (lower complexity AND lower error than everything before it), full list
(10 points, none dominated) -- this IS the complete frontier down to the best error found
in the searched space (8.84e-10 at comp 33; search did not reach 1e-12 within
den<=12/|num|<=48/comp<=40):

```
comp= 5  relerr=3.094e-03  0 + 1*phi^2 + 2*pi*phi
comp= 6  relerr=1.006e-03  1 + 1*pi^2 + 1*pi/phi
comp= 8  relerr=2.455e-04  1/3 + 1*pi^2 + 1*phi^2
comp= 9  relerr=3.229e-05  2 + 1*pi^2 + 3*1/pi
comp=13  relerr=2.604e-05  1/2 + 4/3*phi + 2*pi*phi
comp=14  relerr=9.019e-07  3 + 1*pi^2 + -1/7*1/pi
comp=19  relerr=2.234e-07  -1/6 + 4*pi + 4/3*1/pi
comp=23  relerr=9.187e-08  -2/9 + 7/2*phi^2 + 2*pi/phi
comp=24  relerr=1.108e-08  10 + -4/3*log(phi) + 5*log(2)
comp=33  relerr=8.841e-10  -3 + 21/4*phi^2 + 3*log(2)
```

Candidate C = 15/8 + 4*pi - phi (comp 30) relerr = 6.109e-05 vs D1.
Position on frontier: C is NOT on the Pareto frontier -- it is dominated. In particular
comp=9 (2 + pi^2 + 3/pi, relerr 3.229e-05) is both simpler AND more accurate, and comp=14
(relerr 9.0e-7) is far more accurate at less than half the complexity.

## Test C (delta ~ q * pi^a * phi^b, a in [-2,2], b in [-30,10], q den<=12 |num|<=48)

Top 20 hits by relative error (full list also in testC_hits.txt):
```
a= 2 b=-22 q=22/7   comp(q)=29  relerr=1.605e-04
a=-1 b=-19 q=23     comp(q)=24  relerr=3.527e-04
a=-1 b=-16 q=38/7   comp(q)=45  relerr=5.354e-04
a= 2 b=-20 q=6/5    comp(q)=11  relerr=5.475e-04
a= 0 b=-22 q=31     comp(q)=32  relerr=7.651e-04
a= 1 b=-19 q=7/3    comp(q)=10  relerr=9.110e-04
a= 0 b=-17 q=14/5   comp(q)=19  relerr=9.264e-04
a= 0 b=-16 q=19/11  comp(q)=30  relerr=9.375e-04
a=-2 b= -8 q=4/11   comp(q)=15  relerr=1.153e-03
a= 0 b=-19 q=22/3   comp(q)=25  relerr=1.314e-03
a=-1 b=-17 q=44/5   comp(q)=49  relerr=1.329e-03
a= 2 b=-25 q=40/3   comp(q)=43  relerr=1.340e-03
a= 1 b=-23 q=16     comp(q)=17  relerr=1.355e-03
a= 1 b=-20 q=34/9   comp(q)=43  relerr=1.538e-03
a= 1 b=-17 q=8/9    comp(q)=17  relerr=1.745e-03
a=-1 b=-15 q=37/11  comp(q)=48  relerr=2.025e-03
a= 2 b=-26 q=43/2   comp(q)=45  relerr=2.085e-03
a= 2 b=-24 q=33/4   comp(q)=37  relerr=2.500e-03
a= 1 b=-18 q=13/9   comp(q)=22  relerr=2.553e-03
a=-2 b=-14 q=13/2   comp(q)=15  relerr=2.712e-03
```

No hit anywhere in the search has relerr < 1e-4 AND comp(q) <= 10 -- flagged tier is EMPTY.
(Best relerr overall is 1.605e-04 at comp(q)=29, i.e. high complexity, not a genuine lead.)

Top 10 by complexity(q) among relerr<1e-2 (in testC_hits.txt):
```
a=-2 b=-13 q=4     comp=5  relerr=6.987e-03
a= 0 b=-12 q=1/4   comp=5  relerr=8.889e-03
a= 0 b=-14 q=2/3   comp=5  relerr=9.522e-03
a=-2 b=-12 q=5/2   comp=7  relerr=4.205e-03
a= 2 b=-17 q=2/7   comp=9  relerr=8.035e-03
a=-1 b=-12 q=4/5   comp=9  relerr=9.537e-03
a= 1 b=-19 q=7/3   comp=10 relerr=9.110e-04
a= 2 b=-20 q=6/5   comp=11 relerr=5.475e-04
a= 0 b=-18 q=9/2   comp=11 relerr=5.811e-03
a=-1 b=-10 q=3/10  comp=13 relerr=8.874e-03
```

Specific framework-constant ratios (delta = 7.8336834580104735...e-4):
```
1/(9*phi**6)          = 6.19201000009e-3   delta/this = 0.12651277  relerr(vs delta) = 6.90e+00
3/(2*pi)              = 4.77464829276e-1   delta/this = 0.00164068  relerr(vs delta) = 6.09e+02
1/360                 = 2.77777777778e-3   delta/this = 0.28201260  relerr(vs delta) = 2.55e+00
1/phi**15             = 7.33137435857e-4   delta/this = 1.06851500  relerr(vs delta) = 6.41e-02
1/(4*pi*phi**6)       = 4.43470049635e-3   delta/this = 0.17664515  relerr(vs delta) = 4.66e+00
(15/8-phi)/phi**12    = 7.98038787749e-4   delta/this = 0.98161688  relerr(vs delta) = 1.87e-02
```
None match delta closely; closest are 1/phi**15 (6.4% off) and (15/8-phi)/phi**12 (1.9% off).

## Test D (PSLQ calibration, mpmath.pslq, maxcoeff=1e8, maxsteps=1e7)

```
dps=40  vec=[D1,1,pi,phi]                                   relation=[-5,-3392,0,2136]  height=3392
dps=40  vec=[D1,1,pi,phi,pi*phi,pi^2,phi^2]                  relation=[0,1,0,1,0,0,-1]   height=1
dps=80  vec=[D1,1,pi,phi]              (tol loosened to 1e-30) relation=[-5,-3392,0,2136]  height=3392
dps=80  vec=[D1,1,pi,phi,pi*phi,pi^2,phi^2] (tol loosened to 1e-30) relation=[0,1,0,1,0,0,-1] height=1
```
Note: at dps=80 with default tol, pslq returned None for vec1 (default tolerance too tight
relative to dps); explicitly loosening tol to 1e-30 recovered the same relation as dps=40.

Both relations found are NOT "fake" pi-relations:
- vec1's relation -5*D1-3392+2136*phi=0 is the TRUE algebraic identity already given
  (D1=(2136phi-3392)/5) -- pi's coefficient is 0. This is pslq correctly re-deriving the
  known identity, not a precision artifact.
- vec2's relation is 1+phi-phi^2=0, the trivial golden-ratio identity (comp/height 1),
  entirely independent of D1 and pi -- an artifact of adding phi^2 to a basis that already
  contains phi, not evidence about D1 vs pi.

No spurious relation entangling pi with D1 (i.e., with nonzero pi coefficient) was found at
either precision within maxcoeff=1e8. This calibrates as: at this coefficient-height budget,
PSLQ found zero false D1-pi relations, only the two true algebraic ones the vectors trivially
contain.
