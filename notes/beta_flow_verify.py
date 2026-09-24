"""
Verification script for notes/beta_flow_derivation.md.

Independently reconstructs, from the raw Battye-Sutcliffe profile and
Paper I's own thin-torus integral definitions, the self-consistent
feedback coupling b* = beta*/C^2 and the reduced energy functional
across the beta-flow (beta=0 -> beta=beta*) within the fixed Q_H=2
sector. Checks every intermediate quantity against values already
established/published in the framework, then reports the magnitude
gap against the Delta_1 target this was checking for a match to.

No external dependencies (avoids the broken local scipy/numpy ABI);
uses a plain tan-substitution quadrature good to ~1e-9.
"""

import math

PHI = (1 + 5 ** 0.5) / 2


def quad_0_to_inf(f, n=300_000):
    """Integrate f(x) dx over [0, inf) via x = tan(theta)."""
    total = 0.0
    dtheta = (math.pi / 2) / n
    for i in range(n):
        theta = (i + 0.5) * dtheta
        x = math.tan(theta)
        sec2 = 1 + x * x
        total += f(x) * sec2 * dtheta
    return total


def I_fb_numeric(b, n=300_000):
    """I_fb(b) = integral of [kern/(1+b*kern)] * x dx, kern=8/(1+x^2)^2."""
    if b == 0:
        return 4.0
    return quad_0_to_inf(
        lambda x: (8 / (1 + x ** 2) ** 2) / (1 + b * 8 / (1 + x ** 2) ** 2) * x, n
    )


def I_fb_closed(b):
    """Closed form: I_fb(b) = (4/a) * arctan(a), a = sqrt(8b)."""
    if b == 0:
        return 4.0
    a = math.sqrt(8 * b)
    return (4.0 / a) * math.atan(a)


def I_2a_numeric(n=300_000):
    """I_2a = 128 * integral of x^5/(1+x^2)^6 dx (established closed form: 32/15)."""
    return 128 * quad_0_to_inf(lambda x: x ** 5 / (1 + x ** 2) ** 6, n)


def I_4_numeric(n=300_000):
    """I_4 = 64 * integral of x^3/(1+x^2)^6 dx (established closed form: 8/5)."""
    return 64 * quad_0_to_inf(lambda x: x ** 3 / (1 + x ** 2) ** 6, n)


def solve_bstar(I_2a, tol=1e-14):
    """Solve I_fb(b*) = phi * I_2a by bisection."""
    target = PHI * I_2a
    lo, hi = 0.0, 2.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if I_fb_closed(mid) > target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2


def main():
    print("=== Step 1: I_fb(b) closed form vs numeric quadrature ===")
    for b in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]:
        num = I_fb_numeric(b)
        closed = I_fb_closed(b)
        print(f"  b={b:<8} numeric={num:.10f}  closed={closed:.10f}  diff={abs(num-closed):.2e}")

    print("\n=== Step 2: beta-independent integrals vs established values ===")
    I_2a = I_2a_numeric()
    I_4 = I_4_numeric()
    print(f"  I_2a (numeric)  = {I_2a:.10f}   established 32/15 = {32/15:.10f}")
    print(f"  I_4  (numeric)  = {I_4:.10f}   established 8/5  = {8/5:.10f}")
    print(f"  I_4/I_2a        = {I_4/I_2a:.10f}   established BPS ratio 3/4 = {0.75:.10f}")

    print("\n=== Step 3: solve for b* = beta*/C^2 ===")
    bstar = solve_bstar(I_2a)
    print(f"  b* (this calc)  = {bstar:.10f}")
    print(f"  b* (framework, P4:prop:bstar) approx 0.06715")
    print(f"  match: {'YES' if abs(bstar-0.06715) < 1e-4 else 'NO'}")

    print("\n=== Step 4: reduced energy E(b) = I_2a + (3-phi)*I_fb(b) ===")
    mu2 = 3 - PHI
    E0 = I_2a + mu2 * I_fb_closed(0)
    Ebstar = I_2a + mu2 * I_fb_closed(bstar)
    print(f"  E(0)   = {E0:.6f}")
    print(f"  E(b*)  = {Ebstar:.6f}")
    print(f"  established check: E(b*) should equal 2*phi*I_2a = {2*PHI*I_2a:.6f}")
    print(f"  match: {'YES' if abs(Ebstar - 2*PHI*I_2a) < 1e-8 else 'NO'}")

    print("\n=== Step 5: magnitude check against Delta_1 target ===")
    delta_E = E0 - Ebstar
    ln_ratio_E = math.log(E0 / Ebstar)
    target = 12.824119993955094 * 2 * math.pi / 3  # = 26.858774...
    print(f"  E(0) - E(b*)         = {delta_E:.6f}")
    print(f"  ln(E(0)/E(b*))       = {ln_ratio_E:.6f}")
    print(f"  target log(mu_UV/mu_IR) = {target:.6f}")
    print(f"  ratio target/delta_E     = {target/delta_E:.2f}")
    print(f"  ratio target/ln_ratio_E  = {target/ln_ratio_E:.2f}")
    print("  --> both candidates are ~2 orders of magnitude too small.")
    print("      Naive bare-energy-difference version of the beta-flow")
    print("      hypothesis does NOT reproduce Delta_1 without an extra,")
    print("      currently unidentified physical normalization/prefactor.")


if __name__ == "__main__":
    main()
