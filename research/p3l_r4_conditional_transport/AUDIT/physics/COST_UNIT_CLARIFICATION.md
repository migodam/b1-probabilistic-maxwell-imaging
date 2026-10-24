# Cost unit reconciliation

The earlier two reuse job receipts record 8,003 **full-grid** complex RHS. They also executed 6,144 reduced complex-system RHS columns, recorded separately in their profiles. The user's cap counts all complex RHS, so the conservative combined cap charge is **14,147**, with full-grid/reduced breakdown preserved. There were 768 physical full-grid LU and 768 reduced LU, distinct counts.

The parent owns the append-only cap reconciliation. Do not edit or replay the closed job receipts. The earlier `PAID_JOBS_FINAL_COST.json:new_RHS=8003` is historical full-grid scope, not the reconciled all-complex total. A RHS-saving claim cannot hide reduced columns.

The geometry and Heun geometry jobs use booked full-grid complex current/receiver-adjoint systems only. Decoder roots, real generalized eigensystems, real coordinate calculations and SLSQP parameter solves are not additional complex current-system RHS. Any future reduced complex current-system solve must reserve its columns before computation.
