# Possibility-volume engine: implementation ready, physics disabled

The executable is src/possibility_volume.py. Separate generation is src/possibility_volume_generate.py. Root approval and matching explicit enablement are required before either performs physics. The draft is config/possibility_volume_protocol_DRAFT.json, approved_by_root=false.

No physics or parent generation has been run for this new namespace. Syntax compilation and static semantic review only were performed; there is no micro numerical test or hidden pilot.

## Executable interfaces

1. Engine draft/plan: draft writes the unapproved draft; plan prints the arithmetic cost envelope without opening physics.
2. Separate generator: possibility_volume_generate.py --config <approved protocol>. It draws32 original Gaussian66 parents and independent complex noise streams. Every truth state, solver attempt, material and generator prediction stays in the private output directory. Public online files contain exactly observation/sigma/parent.
3. Online: possibility_volume.py parent --parent R9V_FRESH_000 --config <approved protocol>, or all for the32 listed parents. It opens only online inputs, starts its anchor from zeros66, and does not import the generator.
4. Seal: possibility_volume.py seal --config <approved protocol>. This requires every parent COMPLETE record with all48 arm/proposal rows (including unavailable proposals); then writes PREDICTIONS_SEALED.json. Truth scoring is a separate downstream process and remains blocked until that seal.

## Numerical semantics

The fixed host is eliminated exactly for arbitrary nonuniform P1 ROI materials. With the already-paid cache C, ER, O, Fh:

L_R = I256 - diag(delta) C; j_R = L_R^-1 diag(delta) ER;
F = Fh + O j_R; E_R = ER + C j_R.

Each actual forward reserves8 RHS256 and one LU256. Each required full66 J reserves32 adjoint RHS256 against the same live LU. The derivative is O L_R^-1 diag(E_R) D_x chi, packed into512x66 real whitened columns. Full field/current are reconstructed algebraically when persisting selected endpoints. No new576-grid audit or fallback is enabled. The existing setup cost264 and historical complete study2184 full576 RHS are identified as already paid and never recharged.

The anchor reproduces R7 direct-probit TRF settings: zeros66, cap64 fresh states with one reserved restoration slot, exact TRF, x_scale=jac, ftol=xtol=1e-9, gtol=1e-7, true MAP objective. All anchor work is separately recorded.

Cloud order: pair0, all four radii .02/.05/.10/.15 with +/-; then pair1 in that same radius order. Thus K16 is total, four proposals per radius. Every radius starts independently from the same anchor. B/C share the same8 material-weak Gaussian weights; random66 has separate paired native streams. The32 fresh prior/noise streams use distinct declared SeedSequence role tuples in namespace20261030_R9_POSSIBILITY_VOLUME_V1.

B uses anchor rank8 weak frame, its58D orthogonal complement plus the selected weak ray (59D allowed family). C refreshes the current rank8 frame and uses the current weak projection ray stepwise. The random arm uses the stronger ordinary full66 effective-material geometry. Only the signed-progress tangent row is imposed before original radial retraction. Both weak spectrum and GN model retain the original prior weight1.

The GN MAP model and primary volume qualification differ deliberately. The primary endpoint is the lowest actual dataPhi among visited legal, radius/progress-qualified states. DeltaPhi<=2 and DeltaPrior<=10 remain secondary fields. Absolute95% noise-ball membership, original decoder legality, radius2%, and finite progress>=.5r define primary qualification. Raw proposal and every paid visited state are retained. A raw proposal stays selected if no visited eligible state lowers true dataPhi.

Each proposal has at most6 forwards,4 J-bearing states (including its raw state), and at most3 accepted corrections. Rejects retain their costs. A correction can use a cached J from the current selected point even after its LU is released; a new J is never computed on a released LU. Every selected endpoint has a paid J; scientific infeasibility/low rank is recorded rather than silently relaxed.

## Ledger and data schema

One shared physical worker lock, one live256LU, BLAS1, float64/complex128. Append-only RESOURCE_LEDGER.jsonl distinguishes reservations, attempted RHS/LU, completed solves, failure events, array-I/O costs, controlled wall and CPU. The hardcaps are400000 RHS256,12000 LU256,7200 controlled seconds,8GB peak memory. No reservation is refunded.

Worst bound:1536 proposals x (6x8 + 4x32) =270336 RHS256;32 anchors x64x40 =81920; generation256. Total352512 RHS256,11296 LU256,0 new full576 RHS. This bound fits the root's caps; it is not a prediction of actual work.

Each parent has:
- anchor/ANCHOR.npz and ANCHOR_STATUS.json.
- ANCHOR_FRAME.npz.
- <arm>/<proposal>/SELECTED.npz, RESULT.json, VISITED.json, VISITED_COORDINATES.npz.
- COMPLETE.json with48 rows, explicit slots0..15, proposal seed/order manifest, shared anchor cost, initial frame wall/CPU and cache-load overhead.
- Every actual trial remains under jobs/<parent>/states with x/material/prediction, reduced current/field and J if actually computed.

Rows contain primary qualification, absolute Phi, legacy DeltaPhi/Prior, actual material RMS, signed progress, raw qualification, selected_path, paid trial/J counts, actual wall/CPU, arm-order/window metadata and GT_used=false. The shared anchor is present in cloud storage but becomes primary-qualified only if legal and in the absolute noise ball. Unavailable slots stay in the16-slot denominator.

## Remaining protocol review point

Allowing GN to use the selected weak ray permits its scalar coefficient to cross zero. Signed material progress does not by itself prove positive native weak amplitude because Qs is native-orthogonal, not material-orthogonal. The engine currently reports anchor weak-ray defect and preserves original displacement gates; if positive amplitude is required, root must freeze that explicit condition before enabling. Otherwise the exact invariant is a weak line, with its sign reported rather than silently called a positive ray.

The statistical adapter is being authored separately against this schema. No posterior sampling or posterior-accuracy requirement is imposed by this engine.
