# High-contrast actual true3D supplement

Run with project numerical Python and one BLAS thread:

1. `generate_high.py`: freezes the entire contract before random geometry,
   builds two independent continuous non-spherical piecewise-constant objects,
   generates actual dense dyadic DDA observations on6³/8³, and stores exact
   volume-mean labels in `sealed_truth/`.
2. `run_high.py run`: twelve public-only solves. `solve --case ... --method ...`
   is also available for a finite external queue. Its only data reads are the
   frozen protocol and `data_public/{case}.npz` observations/sigma/bounds.
3. `run_high.py seal`: seals final predictions and minimum-paid-Phi CPU-prefix
   choices. The sequential `run` command invokes this after all12 finish.
4. `score_high.py`: requires all12-seal, then reads labels offline, computes
   material errors and receiver discrepancy from saved paid physical currents,
   and writes tables, figures and cost ledger. It performs no new Maxwell LU.

`model_adapter.py` subclasses the existing R8 VectorPR3D primitives. It supplies
explicit `[0,8]` material bounds to the caller's inherited optimization body
and records actual inherited matrix-action counts. Existing shared R8 defaults
are preserved. All currents/fields remain full3D cell-major vector arrays with
one isotropic real scalar chi shared across components and sources.

Methods are original CIE_sourcewise and Dual_sourcewise adaptations, plus AP/DR.
They share the original physical objective/prior/initial point and common
max1000-iteration/60CPU-second corrector. Native baseline120 cycles/inner12 and
projection3000 steps/audit50/320CPU-second allowance remain distinct measured
proposal costs. AP/DR use grouped ports. The inherited baseline coordinate
default is the SVD range of projected IFFT, explicitly not literal2019 paper
reproduction or original Fourier-gamma coordinates. No A+/Global extra arms are
silently added to the four-method frozen cohort.

All geometry volumes are analytic rectangular intersection integrals of a
continuous material definition, including overlapping boxes/cavities. The
8³→6³ case uses the first object's exact same continuous geometry and its exact
6³ volume mean. It is a paired discretization stressor, not a third independent
sample. Coarse high-index DDA remains a finite-model challenge: cell phase is
large and no continuum-accuracy or universal-recovery claim is supported.

Each result keeps `prediction.npz`, every paid material and physical-Phi/time/
counter history, raw proposal selection, shared corrector termination, and full
projection audit states where applicable. `COST_LEDGER.json` retains physical
LU/RHS, setup, fixed-H factors/RHS, current PR/matrix products, material calls,
operator actions and archive timings separately. Timers/counters nest; do not
add overlapping records. Module imports before solve clocks are excluded.
