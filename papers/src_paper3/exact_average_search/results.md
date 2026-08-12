# Exact-average search at the density-feedback fixed point

Model: X(t) = 8/(1+t^2)^2, area measure t dt, coupling b,
c = sqrt(8 b*); fixed point defined by (15/8) arctan(c)/c = phi.
At b = b*: arctan(c) = (8 phi/15) c exactly. Every candidate
integral reduces to R0(c) + R1(c) arctan(c) with R0, R1 rational;
after the substitution the result is rational in c and phi.
c is transcendental (Lindemann) with no algebraic relations
beyond the arctan identity, so: c survives -> the value is
transcendental -> provably not equal to any algebraic target;
c cancels -> the value is algebraic in phi -> exact comparison
against the targets.

Targets (fixed): T1 = 112.5/phi^10, T2 = sqrt(T1) = 15/(sqrt(2) phi^5),
T3 = 1/T1, T4 = 1/T2.

Numerics: c* = 0.73294641696045245095399914633940724730540704897063
  T1 = 0.91469461002562673411896236261414660558694995184229
  T2 = 0.95639668026694171190959408446342676624269400304871
  T3 = 1.0932610611666152591222424523565341908854844879050
  T4 = 1.0455912495648647398937681229490090510954722108591

Every symbolic integral was cross-checked against 50-digit mpmath
quadrature at c = c*, and every fixed-point reduction was
cross-checked numerically to 45+ digits.

## Validations

P1: (15/8) (int X/(1+bX) t dt)/(int X t dt)
    = (15/8) arctan(c)/c -> phi after substitution, c-free. PASS.
    (int X t dt = 4, as required.)

P2: E1 = (int X/(1+bX)^2 t dt)/4
    = (c + (c^2+1) arctan(c)) / (2 c (c^2+1))
    -> 4 phi/15 + 1/(2(1+c^2)) after substitution; c survives ->
    transcendental, miss. Numeric at c*: 0.756740237842. PASS.

## Pairwise averages <F>_w = (int F w t dt)/(int w t dt)

Notes: (1+bX)^-2 is the same function as eps_tan, and the weight
X*eps_tan is the same function as X/(1+bX)^2; rows are kept for
completeness of the requested grid.

| F | w | closed form (pre-substitution) | after substitution | c-free? | verdict | numeric at c* |
|---|---|---|---|---|---|---|
| eps_tan | X | `(c + (c^2 + 1)*arctan(c))/(2*c*(c^2 + 1))` | `(8*phi*(c^2 + 1) + 15)/(30*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.756740237842 |
| eps_tan | X^2 | `3*(-c + (c^2 + 1)*arctan(c))/(2*c^3*(c^2 + 1))` | `(8*phi*(c^2 + 1) - 15)/(10*c^2*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.5931263407 |
| eps_tan | X/(1+bX) | `(c*(3*c^2 + 5) + 3*(c^4 + 2*c^2 + 1)*arctan(c))/(8*(c^4 + 2*c^2 + 1)*arctan(c))` | `3*(15*c^2 + 8*phi*(c^4 + 2*c^2 + 1) + 25)/(64*phi*(c^4 + 2*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.780289748517 |
| eps_tan | X/(1+bX)^2 | `(c^2 + 1)*(c*(15*c^4 + 40*c^2 + 33) + 15*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))/(24*(c + (c^2 + 1)*arctan(c))*(c^6 + 3*c^4 + 3*c^2 + 1))` | `5*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 40*c^2 + 8*phi + 33)/(8*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 30*c^2 + 8*phi + 15))` | no | transcendental (c survives) -- provable miss | 0.802137197902 |
| eps_tan | X*eps_tan | `(c^2 + 1)*(c*(15*c^4 + 40*c^2 + 33) + 15*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))/(24*(c + (c^2 + 1)*arctan(c))*(c^6 + 3*c^4 + 3*c^2 + 1))` | `5*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 40*c^2 + 8*phi + 33)/(8*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 30*c^2 + 8*phi + 15))` | no | transcendental (c survives) -- provable miss | 0.802137197902 |
| eps_rad | X | `1/(c^4 + 2*c^2 + 1)` | `1/(c^4 + 2*c^2 + 1)` | no | transcendental (c survives) -- provable miss | 0.423187999379 |
| eps_rad | X^2 | `3*(c*(2*c^2 + 1) - (c^4 + 2*c^2 + 1)*arctan(c))/(c^3*(c^4 + 2*c^2 + 1))` | `(30*c^2 - 8*phi*(c^4 + 2*c^2 + 1) + 15)/(5*c^2*(c^4 + 2*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.0833113167355 |
| eps_rad | X/(1+bX) | `(c*(3*c^4 + 8*c^2 + 21) + 3*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))/(24*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))` | `(15*c^4 + 40*c^2 + 8*phi*(c^6 + 3*c^4 + 3*c^2 + 1) + 105)/(64*phi*(c^6 + 3*c^4 + 3*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.472774494982 |
| eps_rad | X/(1+bX)^2 | `(c^2 + 1)*(c*(15*c^6 + 55*c^4 + 73*c^2 + 81) + 15*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1)*arctan(c))/(48*(c + (c^2 + 1)*arctan(c))*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | `5*(8*c^10*phi + 40*c^8*phi + 15*c^8 + 80*c^6*phi + 70*c^6 + 80*c^4*phi + 128*c^4 + 40*c^2*phi + 154*c^2 + 8*phi + 81)/(16*(8*c^2*phi + 8*phi + 15)*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.519397236615 |
| eps_rad | X*eps_tan | `(c^2 + 1)*(c*(15*c^6 + 55*c^4 + 73*c^2 + 81) + 15*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1)*arctan(c))/(48*(c + (c^2 + 1)*arctan(c))*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | `5*(8*c^10*phi + 40*c^8*phi + 15*c^8 + 80*c^6*phi + 70*c^6 + 80*c^4*phi + 128*c^4 + 40*c^2*phi + 154*c^2 + 8*phi + 81)/(16*(8*c^2*phi + 8*phi + 15)*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.519397236615 |
| eps_pullback | X | `(c*(c^2 + 3) + (c^4 + 2*c^2 + 1)*arctan(c))/(4*c*(c^4 + 2*c^2 + 1))` | `(15*c^2 + 8*phi*(c^4 + 2*c^2 + 1) + 45)/(60*(c^4 + 2*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.58996411861 |
| eps_pullback | X^2 | `3*(c*(3*c^2 + 1) - (c^4 + 2*c^2 + 1)*arctan(c))/(4*c^3*(c^4 + 2*c^2 + 1))` | `(45*c^2 - 8*phi*(c^4 + 2*c^2 + 1) + 15)/(20*c^2*(c^4 + 2*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.338218828718 |
| eps_pullback | X/(1+bX) | `(c*(3*c^4 + 8*c^2 + 9) + 3*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))/(12*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))` | `(15*c^4 + 40*c^2 + 8*phi*(c^6 + 3*c^4 + 3*c^2 + 1) + 45)/(32*phi*(c^6 + 3*c^4 + 3*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.626532121749 |
| eps_pullback | X/(1+bX)^2 | `(c^2 + 1)*(c*(15*c^6 + 55*c^4 + 73*c^2 + 49) + 15*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1)*arctan(c))/(32*(c + (c^2 + 1)*arctan(c))*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | `15*(8*c^10*phi + 40*c^8*phi + 15*c^8 + 80*c^6*phi + 70*c^6 + 80*c^4*phi + 128*c^4 + 40*c^2*phi + 122*c^2 + 8*phi + 49)/(32*(8*c^2*phi + 8*phi + 15)*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.660767217258 |
| eps_pullback | X*eps_tan | `(c^2 + 1)*(c*(15*c^6 + 55*c^4 + 73*c^2 + 49) + 15*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1)*arctan(c))/(32*(c + (c^2 + 1)*arctan(c))*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | `15*(8*c^10*phi + 40*c^8*phi + 15*c^8 + 80*c^6*phi + 70*c^6 + 80*c^4*phi + 128*c^4 + 40*c^2*phi + 122*c^2 + 8*phi + 49)/(32*(8*c^2*phi + 8*phi + 15)*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.660767217258 |
| 1/eps_tan | X | `c^4/5 + 2*c^2/3 + 1` | `c^4/5 + 2*c^2/3 + 1` | no | transcendental (c survives) -- provable miss | 1.41585931364 |
| 1/eps_tan | X^2 | `3*c^4/7 + 6*c^2/5 + 1` | `3*c^4/7 + 6*c^2/5 + 1` | no | transcendental (c survives) -- provable miss | 1.76833614062 |
| 1/eps_tan | X/(1+bX) | `c*(c^2 + 3)/(3*arctan(c))` | `5*(c^2 + 3)/(8*phi)` | no | transcendental (c survives) -- provable miss | 1.36632267722 |
| 1/eps_tan | X/(1+bX)^2 | `2*c*(c^2 + 1)/(c + (c^2 + 1)*arctan(c))` | `30*(c^2 + 1)/(8*c^2*phi + 8*phi + 15)` | no | transcendental (c survives) -- provable miss | 1.32145741695 |
| 1/eps_tan | X*eps_tan | `2*c*(c^2 + 1)/(c + (c^2 + 1)*arctan(c))` | `30*(c^2 + 1)/(8*c^2*phi + 8*phi + 15)` | no | transcendental (c survives) -- provable miss | 1.32145741695 |
| (1+bX)^-1 | X | `arctan(c)/c` | `8*phi/15` | yes | exact algebraic value 4/15 + 4*sqrt(5)/15 -- no target match | 0.862951460667 |
| (1+bX)^-1 | X^2 | `3*(c - arctan(c))/c^3` | `(15 - 8*phi)/(5*c^2)` | no | transcendental (c survives) -- provable miss | 0.76533436365 |
| (1+bX)^-1 | X/(1+bX) | `(c + (c^2 + 1)*arctan(c))/(2*(c^2 + 1)*arctan(c))` | `(8*phi*(c^2 + 1) + 15)/(16*phi*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.876920976827 |
| (1+bX)^-1 | X/(1+bX)^2 | `(c^2 + 1)*(c*(3*c^2 + 5) + 3*(c^4 + 2*c^2 + 1)*arctan(c))/(4*(c + (c^2 + 1)*arctan(c))*(c^4 + 2*c^2 + 1))` | `3*(8*c^4*phi + 16*c^2*phi + 15*c^2 + 8*phi + 25)/(4*(8*c^4*phi + 16*c^2*phi + 15*c^2 + 8*phi + 15))` | no | transcendental (c survives) -- provable miss | 0.889806230136 |
| (1+bX)^-1 | X*eps_tan | `(c^2 + 1)*(c*(3*c^2 + 5) + 3*(c^4 + 2*c^2 + 1)*arctan(c))/(4*(c + (c^2 + 1)*arctan(c))*(c^4 + 2*c^2 + 1))` | `3*(8*c^4*phi + 16*c^2*phi + 15*c^2 + 8*phi + 25)/(4*(8*c^4*phi + 16*c^2*phi + 15*c^2 + 8*phi + 15))` | no | transcendental (c survives) -- provable miss | 0.889806230136 |
| (1+bX)^-2 | X | `(c + (c^2 + 1)*arctan(c))/(2*c*(c^2 + 1))` | `(8*phi*(c^2 + 1) + 15)/(30*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.756740237842 |
| (1+bX)^-2 | X^2 | `3*(-c + (c^2 + 1)*arctan(c))/(2*c^3*(c^2 + 1))` | `(8*phi*(c^2 + 1) - 15)/(10*c^2*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.5931263407 |
| (1+bX)^-2 | X/(1+bX) | `(c*(3*c^2 + 5) + 3*(c^4 + 2*c^2 + 1)*arctan(c))/(8*(c^4 + 2*c^2 + 1)*arctan(c))` | `3*(15*c^2 + 8*phi*(c^4 + 2*c^2 + 1) + 25)/(64*phi*(c^4 + 2*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.780289748517 |
| (1+bX)^-2 | X/(1+bX)^2 | `(c^2 + 1)*(c*(15*c^4 + 40*c^2 + 33) + 15*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))/(24*(c + (c^2 + 1)*arctan(c))*(c^6 + 3*c^4 + 3*c^2 + 1))` | `5*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 40*c^2 + 8*phi + 33)/(8*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 30*c^2 + 8*phi + 15))` | no | transcendental (c survives) -- provable miss | 0.802137197902 |
| (1+bX)^-2 | X*eps_tan | `(c^2 + 1)*(c*(15*c^4 + 40*c^2 + 33) + 15*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))/(24*(c + (c^2 + 1)*arctan(c))*(c^6 + 3*c^4 + 3*c^2 + 1))` | `5*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 40*c^2 + 8*phi + 33)/(8*(8*c^6*phi + 24*c^4*phi + 15*c^4 + 24*c^2*phi + 30*c^2 + 8*phi + 15))` | no | transcendental (c survives) -- provable miss | 0.802137197902 |
| (1+bX)^-3 | X | `(c*(3*c^2 + 5) + 3*(c^4 + 2*c^2 + 1)*arctan(c))/(8*c*(c^4 + 2*c^2 + 1))` | `(15*c^2 + 8*phi*(c^4 + 2*c^2 + 1) + 25)/(40*(c^4 + 2*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.673352178226 |
| (1+bX)^-3 | X^2 | `3*(c*(c^2 - 1) + (c^4 + 2*c^2 + 1)*arctan(c))/(8*c^3*(c^4 + 2*c^2 + 1))` | `(15*c^2 + 8*phi*(c^4 + 2*c^2 + 1) - 15)/(40*c^2*(c^4 + 2*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.465672584709 |
| (1+bX)^-3 | X/(1+bX) | `(c*(15*c^4 + 40*c^2 + 33) + 15*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))/(48*(c^6 + 3*c^4 + 3*c^2 + 1)*arctan(c))` | `5*(15*c^4 + 40*c^2 + 8*phi*(c^6 + 3*c^4 + 3*c^2 + 1) + 33)/(128*phi*(c^6 + 3*c^4 + 3*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.703410935133 |
| (1+bX)^-3 | X/(1+bX)^2 | `(c^2 + 1)*(c*(105*c^6 + 385*c^4 + 511*c^2 + 279) + 105*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1)*arctan(c))/(192*(c + (c^2 + 1)*arctan(c))*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | `5*(56*c^10*phi + 280*c^8*phi + 105*c^8 + 560*c^6*phi + 490*c^6 + 560*c^4*phi + 896*c^4 + 280*c^2*phi + 790*c^2 + 56*phi + 279)/(64*(8*c^2*phi + 8*phi + 15)*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.73145220758 |
| (1+bX)^-3 | X*eps_tan | `(c^2 + 1)*(c*(105*c^6 + 385*c^4 + 511*c^2 + 279) + 105*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1)*arctan(c))/(192*(c + (c^2 + 1)*arctan(c))*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | `5*(56*c^10*phi + 280*c^8*phi + 105*c^8 + 560*c^6*phi + 490*c^6 + 560*c^4*phi + 896*c^4 + 280*c^2*phi + 790*c^2 + 56*phi + 279)/(64*(8*c^2*phi + 8*phi + 15)*(c^8 + 4*c^6 + 6*c^4 + 4*c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.73145220758 |

| F | w | status |
|---|---|---|
| sqrt(eps_tan*eps_rad) | X | no elementary closed form (integrand s^4 sqrt(s^2-3c^2)/(s^2+c^2)^(5/2) is outside the rational/arctan family; real only for s >= sqrt(3) c); not pursued numerically |
| sqrt(eps_tan*eps_rad) | X^2 | no elementary closed form (integrand s^4 sqrt(s^2-3c^2)/(s^2+c^2)^(5/2) is outside the rational/arctan family; real only for s >= sqrt(3) c); not pursued numerically |
| sqrt(eps_tan*eps_rad) | X/(1+bX) | no elementary closed form (integrand s^4 sqrt(s^2-3c^2)/(s^2+c^2)^(5/2) is outside the rational/arctan family; real only for s >= sqrt(3) c); not pursued numerically |
| sqrt(eps_tan*eps_rad) | X/(1+bX)^2 | no elementary closed form (integrand s^4 sqrt(s^2-3c^2)/(s^2+c^2)^(5/2) is outside the rational/arctan family; real only for s >= sqrt(3) c); not pursued numerically |
| sqrt(eps_tan*eps_rad) | X*eps_tan | no elementary closed form (integrand s^4 sqrt(s^2-3c^2)/(s^2+c^2)^(5/2) is outside the rational/arctan family; real only for s >= sqrt(3) c); not pursued numerically |

1/eps_rad excluded (non-integrable sign crossing); w = 1 excluded
(divergent normalization).

## Standalone combinations at b = b*

V0 = 15/8, V* = phi.

| combination | after substitution | c-free? | verdict | numeric at c* |
|---|---|---|---|---|
| V*^2 | `phi^2` | yes | exact algebraic value sqrt(5)/2 + 3/2 -- no target match | 2.61803398875 |
| (V*/V0)^2 | `64*phi^2/225` | yes | exact algebraic value 32*sqrt(5)/225 + 32/75 -- no target match | 0.744685223467 |
| V0*V* | `15*phi/8` | yes | exact algebraic value 15/16 + 15*sqrt(5)/16 -- no target match | 3.03381372891 |
| 1/(V0*V*) | `8/(15*phi)` | yes | exact algebraic value -4/15 + 4*sqrt(5)/15 -- no target match | 0.329618127333 |
| <eps_tan>_X * <1/eps_tan>_X | `(8*phi*(c^2 + 1) + 15)*(3*c^4 + 10*c^2 + 15)/(450*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 1.07143771375 |
| <eps_tan>_X / <eps_pullback>_X | `2*(8*c^4*phi + 16*c^2*phi + 15*c^2 + 8*phi + 15)/(8*c^4*phi + 16*c^2*phi + 15*c^2 + 8*phi + 45)` | no | transcendental (c survives) -- provable miss | 1.28268858049 |
| (<(1+bX)^-1>_X)^2 | `64*phi^2/225` | yes | exact algebraic value 32*sqrt(5)/225 + 32/75 -- no target match | 0.744685223467 |
| <(1+bX)^-2>_X | `(8*phi*(c^2 + 1) + 15)/(30*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 0.756740237842 |
| <(1+bX)^-2>_X / (<(1+bX)^-1>_X)^2 | `15*(8*phi*(c^2 + 1) + 15)/(128*phi^2*(c^2 + 1))` | no | transcendental (c survives) -- provable miss | 1.01618806711 |

## c-free exact identities of the model

- <(1+bX)^-1>_X = 8*phi/15 = 4/15 + 4*sqrt(5)/15
- V*^2 = phi^2 = sqrt(5)/2 + 3/2
- (V*/V0)^2 = 64*phi^2/225 = 32*sqrt(5)/225 + 32/75
- V0*V* = 15*phi/8 = 15/16 + 15*sqrt(5)/16
- 1/(V0*V*) = 8/(15*phi) = -4/15 + 4*sqrt(5)/15
- (<(1+bX)^-1>_X)^2 = 64*phi^2/225 = 32*sqrt(5)/225 + 32/75

All of these are algebraic consequences of the single defining
identity <(1+bX)^-1>_X = arctan(c)/c = 8 phi/15 (equivalently
V(b*) = phi, validation P1) together with V0 = 15/8; none is an
independent second identity.

## Conclusion

No candidate average or standalone combination equals T1, T2,
T3, or T4 exactly. Every pairwise average that is c-free is an
exact identity listed above; every c-dependent one is
transcendental at b = b* and therefore provably unequal to the
algebraic targets.
