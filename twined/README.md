# Twined indices in the Z2 and Z3 CHL models (Part II)

Code and logs for Part II of the series *Positivity of Black Hole Indices in CHL Models*: the sector indices
S_j = (1/N) sum_k w^{-jk} d_{*,g^k} of single-centred black holes in the Z_N x Z_N twisted-twining models, N = 2, 3,
are positive for every charge, and |d_{*,g}| < 0.03 d_*.

The untwined values d_* come from Part I (folder `allcharges/`). `allcharges/code/sweep_3_2500_0.json` holds the exact d_* for N = 3
and every reduced charge with D <= 2500 (73,147 charges, all positive; log in `allcharges/logs/blocks_3_2500.log`).

## Code (`code/`)

| file | what it does |
|---|---|
| `jac.py`, `cD.py` | exact coefficients c_N(Delta) of 2 phi_{0,2} and 2 phi_{0,3} (pickles `cD_2_400.pkl`, `cD_3_300.pkl`) |
| `tt.py`, `exact_tw.py` | exact expansion of 1/Phi_N^g in the chamber C_0 (integer arithmetic) |
| `smallcheck.py` | twined index for every reduced charge with D <= 400, against exact d_* |
| `lemmaA2.py` | the constant A_N of the exponent bound, with the tail bound |
| `lam.py`, `rsup.py`, `rsup2.py` | the remainder R(Y) and its uniform bound over the 82 cells of the family Y(r) |
| `tailG.py` | the bound U_N(m,n,l;r) and the lower bound L_N(m,D) from Part I |
| `midcheck.py` | N = 3, 400 < D <= 2500: U_3 against exact d_* |
| `boxcheck.py`, `boxcheck2.py` | box certification for D_0 < D <= 1e8 (boxcheck2 records the certified margin) |
| `tailcheck.py` | constants of the closed-form estimate for D > 1e8 |
| `dualcheck.py` | duality invariance of the twined index on 368 pairs of charges |
| `g2theta.py`, `match.py` | the genus-two theta identity for N = 2 |

Run from `code/`, e.g. `python3 rsup2.py 2 1.3 0.15`, `python3 boxcheck2.py 2 1.3 0.15 1.23 400 1e8 3.6`,
`python3 midcheck.py ../../allcharges/code/sweep_3_2500_0.json`. Requires numpy, scipy and python-flint.

## Logs (`logs/`)

| log | result |
|---|---|
| `lemmaA.log` | A_2 <= 11.40800, A_3 <= 6.66594 (used: 11.5, 6.7) |
| `rsup2_2.log`, `rsup2_3.log` | R(Y(r)) <= 1.2201 (N=2), 2.0699 (N=3) on every cell (used: 1.23, 2.08) |
| `small_2.log`, `small_3.log` | D <= 400: 2,661 and 5,168 reduced charges, max abs(d_g)/d_* = 0.0222 and 0.0294 |
| `mid_3.log` | N = 3, 400 < D <= 2500: 67,979 charges, minimal margin 16.32 |
| `box2_2.log`, `box2_3.log` | D_0 < D <= 1e8: 31,840 and 36,538 boxes, certified margin >= 3.60 |
| `tailcheck.log` | D > 1e8: margins above 4000 |
| `dualcheck.log` | 276 + 92 pairs, no mismatch |
| `theta_identity.log` | theta identity for N = 2 through q^5 s^5 |

`data/tw_2_400.json`, `data/tw_3_400.json`: entries `[m, n, l, D, d_*, d_{*,g}]` for every reduced charge with D <= 400.
