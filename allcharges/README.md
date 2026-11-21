# All charges of the Z3, Z5 and Z7 CHL models

Programs and output logs for the paper "Positive Single-Centred Black Hole Indices for All Charges in the Z3, Z5 and Z7 CHL Models".

| Step | Program | Output |
|---|---|---|
| Exact indices for small discriminants | `fastinit.py`, `exact.py`, `sweep.py`, `runblocks.py` | `logs/blocks_*.log` |
| Gamma_0(N) invariance check | `gamma0.py` | `logs/gamma0_check.log` |
| Wall formula for the polar data | `walls.py`, `wallsI.py` | |
| Exact expansion against exact values | `rad.py`, `radK.py` | `logs/radK*.log` |
| Remaining terms, per charge (m < m_U) | `errterms.py`, `errrun.py` | `logs/err_*.log` |
| Remaining terms, uniform bound | `uniform.py`, `crude.py` | `logs/uniform_*.log` |
| Away from the corner | `termwise_sup.py` | `logs/tsup_*.log` |
| Corner, m <= 400, per charge | `cornerL.py`, `cornertail.py` | `logs/cornerL_*_all.log`, `logs/cornertail_*.log` |
| Corner, m >= 300 | `hexp.py`, `largem_E.py`, `largem_grid.py` | `logs/Ecells_*`, `logs/Gcells_*` |

The wall bounds in `errterms.py` and `uniform.py` use the weight m+2-a-b; the paper multiplies them by sqrt(2), which covers the general wall weight |l_G| <= sqrt(2)(m+2).

Requirements: Python 3, numpy, sympy, mpmath, python-flint (for exact series and Arb ball arithmetic). Environment variable `D0` sets the discriminant bound of the exact computation (400 for N=5,7 and 2500 for N=3).

## Proofs of the structural inputs (appendices A–C)

- `fixline.py`: the fixed-curve lattice search for the elliptic elements of Γ0(5) and Γ0(7) (Appendix A); output in `logs/fixline.log`.
- `gamma0e.py`: direct check that d* is invariant under those elliptic elements (`logs/gamma0e_5.log`, `logs/gamma0e_7.log`).
- `myformula.py`: the corrected exact expansion of Theorem 4.3, compared with exact values (`logs/my*.log`, `logs/myT*.log`).
- `knew.py`, `kn3.py`: bounds for the K-terms of the corrected expansion (`logs/knew_*.log`); `crude.py` is updated to the same form.
