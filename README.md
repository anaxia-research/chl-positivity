# chl-positivity

Code and verification logs for **Positivity of Black Hole Indices in CHL Models I: The Four Prime Orders** by Deep Bhattacharjee.

The paper proves Sen's positivity conjecture for every single-centred quarter-BPS black hole in the Z2, Z3, Z5 and Z7 CHL models by one argument, gives two-sided bounds that fix each index up to a bounded factor, and writes the first row Q^2 = 2/N in closed form through the theta series of a maximal quaternion order.

## Layout

- `allcharges/`: the all-charges proof (exact expansion, error bounds, corner certification, two-sided bounds), with its output logs in `allcharges/logs/`. See `allcharges/README.md`.
- `code/`: the exact first-row computations (numerator identity, quaternion theta identity, closed form) and the shared engine `chl.py`.

## First row (run inside `code/`)

Requirements: Python 3.10 or newer, and `pip install python-flint sympy numpy mpmath matplotlib`.

| Command | What it checks | Log |
|---|---|---|
| `python3 verify_num.py N X` | The numerator identity and the quaternion theta identity through x^X, for N = 2, 3, 5, 7 | `allcharges/logs/num_theta.log` |
| `python3 verify_closed.py N X` | The closed form of g(sigma/N) phi^F_1 against the finite part computed from the product formula | `allcharges/logs/closed_N.log` |
| `python3 verify_W.py X` | The forms W_N, their eta-quotient expressions and Eisenstein-cusp decompositions through q^X | `allcharges/logs/W_forms.log` |
| `python3 proof_row1.py` | The inequalities in the sign argument for N = 3, 5, 7 and the small cases | |
| `python3 box.py N X M` | Exact check of every charge with m <= M and m <= n < X | |

`chl.py` expands the product formula in exact rational arithmetic and computes polar and finite parts; `quat.py`, `lattice.py` build the maximal orders (Hurwitz quaternions for N = 2) and their theta series; `row1_lib.py` computes Hurwitz class numbers.

The paper uses: `verify_num.py N 40`, `verify_closed.py N 62`, `verify_W.py 300`.
