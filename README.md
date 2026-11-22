# chl-positivity

Code for **Positive First-Row Black Hole Indices in CHL Models** by Deep Bhattacharjee.

The paper proves Sen's positivity conjecture for every single-centred black hole with electric charge Q² = 2/N in the Z3, Z5 and Z7 CHL models. It also gives an exact closed formula for the finite part of the first Fourier–Jacobi coefficient of the inverse dyon partition function.

## Layout

- `figures/`: created by `figs.py`.
- `code/`: the programs behind every computation in the paper.

## Running the checks

Requirements: Python 3.10 or newer, and `pip install python-flint sympy numpy mpmath matplotlib`.

Run every command from inside `code/`.

| Command | What it checks |
|---|---|
| `python3 verify_id.py N X` | Proposition 3.1 (numerator identity) through x^X. Use X ≥ 12, 60, 168 for N = 3, 5, 7. |
| `python3 match2.py N X` | Proposition 4.1 (quaternion theta identity) through x^X. |
| `python3 verify_K.py N X` | Theorem 1.2: the closed form against the finite part computed from the product formula. |
| `python3 W.py` | The forms W_N and their Eisenstein–cusp decompositions. |
| `python3 proof_row1.py` | The inequalities in the proof of Proposition 6.2 and the list of small cases that are read off exactly. |
| `python3 box.py N X M` | The exact check of every charge with m ≤ M and m ≤ n < X (Section 7). |
| `python3 figs.py` | Redraws the figures of the paper into `figures/`. |

`chl.py` is the core engine. It expands the product formula in exact rational arithmetic and computes the polar and finite parts. `quat.py`, `lattice.py` and `lat2.py` build the maximal quaternion orders and their theta series. `row1_lib.py` computes Hurwitz class numbers.

The paper uses these settings:
- the identities through x^40 (N = 3), x^70 (N = 5) and x^175 (N = 7);
- the theta identity through x^14, x^62 and x^170;
- the box check `box.py N 121 6`.
