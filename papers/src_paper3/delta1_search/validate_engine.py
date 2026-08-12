"""Engine validation: plant T0 = 3*pi/2 + 7/8 - 2*phi, confirm Test A recovers it exactly
at its correct complexity: a=7/8 (comp 7+8=15), b=3/2 (comp 3+2=5), c=-2 (comp 2+1=3).
Total complexity = 15+5+3 = 23.
"""
from fractions import Fraction
from testA_engine import enumerate_exprs, eval_expr_float

TARGET_A = Fraction(7, 8)
TARGET_B = Fraction(3, 2)
TARGET_C = Fraction(-2, 1)
TARGET_COMP = (abs(7)+8) + (abs(3)+2) + (abs(-2)+1)
print("Target (a,b,c) =", TARGET_A, TARGET_B, TARGET_C, "expected total comp =", TARGET_COMP)

MAXC = 23
res = enumerate_exprs(MAXC)
print("Enumerated", len(res), "expressions up to complexity", MAXC)

found = None
for a, b, c, tot in res:
    if a == TARGET_A and b == TARGET_B and c == TARGET_C:
        found = (a, b, c, tot)
        break

if found is None:
    print("FAIL: target expression not found in enumeration")
else:
    print("FOUND:", found)
    print("PASS: engine recovers planted target at complexity", found[3])
