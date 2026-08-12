# SU(2)_3 modular data: exact verification and T_UV naturalness test
q = exp(i*pi/5). Labels j in {0, 1/2, 1, 3/2}. All checks are exact
sympy simplifications to 0, cross-checked at 50-digit numeric precision.

## Quantum dimensions and D^2

| quantity | exact value | numeric (30 dp) |
|---|---|---|
| d_0 | 1 | 1.0 |
| d_1/2 | 1/2 + sqrt(5)/2 | 1.61803398874989484820458683437 |
| d_1 | 1/2 + sqrt(5)/2 | 1.61803398874989484820458683437 |
| d_3/2 | 1 | 1.0 |
| D^2 | sqrt(5) + 5 | 7.23606797749978969640917366873 |
| D^2 = 2*phi*sqrt(5)? | -- | MATCH |

## Twists h_j = j(j+1)/5

| j | h_j |
|---|---|
| 0 | 0 |
| 1/2 | 3/20 |
| 1 | 2/5 |
| 3/2 | 3/4 |

## S-matrix (S_ab = sqrt(2/5) sin((2a+1)(2b+1) pi/5))

| a\\b | 0 | 1/2 | 1 | 3/2 |
|---|---|---|---|---|
| 0 | sqrt(1/4 - sqrt(5)/20) | sqrt(sqrt(5)/20 + 1/4) | sqrt(sqrt(5)/20 + 1/4) | sqrt(1/4 - sqrt(5)/20) |
| 1/2 | sqrt(sqrt(5)/20 + 1/4) | sqrt(1/4 - sqrt(5)/20) | -sqrt(1/4 - sqrt(5)/20) | -sqrt(sqrt(5)/20 + 1/4) |
| 1 | sqrt(sqrt(5)/20 + 1/4) | -sqrt(1/4 - sqrt(5)/20) | -sqrt(1/4 - sqrt(5)/20) | sqrt(sqrt(5)/20 + 1/4) |
| 3/2 | sqrt(1/4 - sqrt(5)/20) | -sqrt(sqrt(5)/20 + 1/4) | sqrt(sqrt(5)/20 + 1/4) | -sqrt(1/4 - sqrt(5)/20) |

## Monodromy of 1/2 x 1/2 -> {0, 1}

| channel c | angle (deg) | angle mod 360 |
|---|---|---|
| 0 | -108 | 252 |
| 1 | 36 | 36 |

Order of q = exp(i*pi/5): 10

## L1 and L2

- L1: 10 * 36 = 360 (order of q, Q = 10); L1 holds: YES
- L2: p_0 = 3/2 - sqrt(5)/2, p_1 = -1/2 + sqrt(5)/2, p_0+p_1 = 1
- L2: 360*p_0 = 540 - 180*sqrt(5) = golden angle (numeric 137.507764050037854646348739628 deg)

## T_UV identities

- T_UV = 16*phi^8/5 = phi^6/sin(pi/5)^4 exactly: YES
- T_UV numeric (30 dp): 150.331884043992933799348235269
- S_00 = sqrt(2/5)*sin(pi/5), exact: sqrt(1/4 - sqrt(5)/20)
- Unique S_00-form found by scan: T_UV = 4/25 * phi^6 / S_00^4
- Verified identity: T_UV = (4/25) * phi^6 / S_00^4 : YES

## Pre-registered vocabulary test against T_UV

Vocabulary: {S_ab, 1/S_ab, d_j, D, D^2, theta_j, p_c, sin(pi/5)}. Only the forms explicitly listed below were tested; no free search.

| form | equals T_UV exactly? |
|---|---|
| 1/S_00^4 | no |
| phi^0/S_00^2 | no |
| phi^1/S_00^2 | no |
| phi^2/S_00^2 | no |
| phi^3/S_00^2 | no |
| phi^4/S_00^2 | no |
| phi^5/S_00^2 | no |
| phi^6/S_00^2 | no |
| phi^7/S_00^2 | no |
| phi^8/S_00^2 | no |
| (4/25)*phi^6/S_00^4  [S_00-form of T_UV] | YES |
| D^2*phi^0/5 | no |
| D^2*phi^1/5 | no |
| D^2*phi^2/5 | no |
| D^2*phi^3/5 | no |
| D^2*phi^4/5 | no |
| D^2*phi^5/5 | no |
| D^2*phi^6/5 | no |
| D^2*phi^7/5 | no |
| D^2*phi^8/5 | no |
| (phi/S_00)^2 | no |
| (phi/S_00)^4/D^2 | no |
| p_0^-0 * p_1^-0 | no |
| p_0^-0 * p_1^-1 | no |
| p_0^-0 * p_1^-2 | no |
| p_0^-0 * p_1^-3 | no |
| p_0^-0 * p_1^-4 | no |
| p_0^-0 * p_1^-5 | no |
| p_0^-0 * p_1^-6 | no |
| p_0^-1 * p_1^-0 | no |
| p_0^-1 * p_1^-1 | no |
| p_0^-1 * p_1^-2 | no |
| p_0^-1 * p_1^-3 | no |
| p_0^-1 * p_1^-4 | no |
| p_0^-1 * p_1^-5 | no |
| p_0^-1 * p_1^-6 | no |
| p_0^-2 * p_1^-0 | no |
| p_0^-2 * p_1^-1 | no |
| p_0^-2 * p_1^-2 | no |
| p_0^-2 * p_1^-3 | no |
| p_0^-2 * p_1^-4 | no |
| p_0^-2 * p_1^-5 | no |
| p_0^-2 * p_1^-6 | no |
| p_0^-3 * p_1^-0 | no |
| p_0^-3 * p_1^-1 | no |
| p_0^-3 * p_1^-2 | no |
| p_0^-3 * p_1^-3 | no |
| p_0^-3 * p_1^-4 | no |
| p_0^-3 * p_1^-5 | no |
| p_0^-3 * p_1^-6 | no |
| p_0^-4 * p_1^-0 | no |
| p_0^-4 * p_1^-1 | no |
| p_0^-4 * p_1^-2 | no |
| p_0^-4 * p_1^-3 | no |
| p_0^-4 * p_1^-4 | no |
| p_0^-4 * p_1^-5 | no |
| p_0^-4 * p_1^-6 | no |
| p_0^-5 * p_1^-0 | no |
| p_0^-5 * p_1^-1 | no |
| p_0^-5 * p_1^-2 | no |
| p_0^-5 * p_1^-3 | no |
| p_0^-5 * p_1^-4 | no |
| p_0^-5 * p_1^-5 | no |
| p_0^-5 * p_1^-6 | no |
| p_0^-6 * p_1^-0 | no |
| p_0^-6 * p_1^-1 | no |
| p_0^-6 * p_1^-2 | no |
| p_0^-6 * p_1^-3 | no |
| p_0^-6 * p_1^-4 | no |
| p_0^-6 * p_1^-5 | no |
| p_0^-6 * p_1^-6 | no |
| 360*p_0/p_1^1 | no |
| 360*p_0/p_1^2 | no |
| 360*p_0/p_1^3 | no |
| 360*p_0/p_1^4 | no |

**Exact hits among tested vocabulary forms: 1** -> (4/25)*phi^6/S_00^4  [S_00-form of T_UV]

Note: the S_00-form identity T_UV = (4/25)*phi^6/S_00^4 (Section 'T_UV identities' above) is an exact restatement of the definition of T_UV in the allowed vocabulary, listed separately from the vocabulary scan since it is a substitution identity (sin(pi/5) -> S_00), not a coincidental match of an independent form.

## T_UV / T_attractor identities

- T_UV/T_attractor = phi^10/112.5, verified exactly: YES
- T_attractor/T_UV = (1/2)*(15/phi^5)^2, verified exactly: YES
- T_attractor identified as 360*p_0 (golden angle from L2), verified exactly: YES
- T_attractor = 360/phi^2, numeric (30 dp): 137.507764050037854646348739628
