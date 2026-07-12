import mpmath as mp
from common import D1

def run_pslq(dps, maxcoeff, vecnames_vals):
    mp.mp.dps = dps
    vec = [f() for name, f in vecnames_vals]
    rel = mp.pslq(vec, maxcoeff=maxcoeff, maxsteps=10**6)
    return rel

specs = [
    ("D1,1,pi,phi", lambda: [D1, mp.mpf(1), mp.pi, (1+mp.sqrt(5))/2]),
    ("D1,1,pi,phi,pi*phi,pi^2,phi^2", lambda: [D1, mp.mpf(1), mp.pi, (1+mp.sqrt(5))/2,
                                                mp.pi*((1+mp.sqrt(5))/2), mp.pi**2, ((1+mp.sqrt(5))/2)**2]),
]

results = []
for dps in (40, 80):
    for maxcoeff in (10**8,):
        for name, f in specs:
            mp.mp.dps = dps
            vec = f()
            rel = mp.pslq(vec, maxcoeff=maxcoeff, maxsteps=10**6)
            height = max(abs(c) for c in rel) if rel else None
            results.append((dps, maxcoeff, name, rel, height))

with open("testD_pslq.txt", "w") as f:
    for dps, maxcoeff, name, rel, height in results:
        f.write(f"dps={dps} maxcoeff={maxcoeff:.0e} vec=[{name}]\n")
        f.write(f"  relation: {rel}\n")
        f.write(f"  coefficient height: {height}\n\n")

for r in results:
    print(r)
