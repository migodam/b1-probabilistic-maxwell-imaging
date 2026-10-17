# F3-first entrypoint: audit, compete, then conditionally authorize physics

RUN_MODE = audit_only
MAX_NEW_FULL_STATES = 0
MAX_NEW_COMPLEX_RHS = 0
MAX_NEW_GPU_TRAINING = 0
MAX_NEW_MIE = 0
MAX_NEW_GENERAL_3D = 0

Read in this order:
1. `integration/F3_A17X_DECISION_20261017.md`
2. `evidence/f3_pinned/F3_EVIDENCE_INDEX.md` and actual pinned F3 files
3. `report/RESEARCH_A17X_FULL_REPORT_ZH.md` (sections 4-7, 9, 11, 13, 15-19)
4. `protocol/CODEX_NEXT_ROUND_ZH.md` (ARCHIVED plan: not blanket run permission)
5. A17 primary theory `CORE_THEORY.tex`, `CONDITIONAL_EXCHANGE_THEORY.tex`.

Run `python checks/check_finite_source_transport.py` and record new CPU-only check provenance. No new physical forward/adjoint/J/GPUs at this stage.

Tasks now:
A. Freeze current F3 facts. Keep former F1 Gate2 FAIL distinct from newly passed F3 F1.5 and completed F3 F2. Do not rerun P0 four arms, original prior MC or older tests; reuse saved endpoints/metrics, label known parent=1.
B. Verify whether P1 is genuinely COMPLETE, has a frozen physically legal prior generator/law, a legitimate imaging task-scale ambiguity witness and reusable Maxwell states. If P1 is pending, report UNKNOWN and prepare protocol without running it.
C. Audit admissible port support U under the ACTUAL P1 prior (or synthetically under P0 if P1 is not ready). A low-dimensional latent is not low-rank DeltaX. For any port compression, quantify operator and receiver truncation errors. Do not silently replace the prior.
D. Check multi-source material consistency, vector block requirements, rational-source prior legality and tail-bound quality by existing cached/synthetic algebra. Mark any real-data stage NOT_RUN if inputs are absent.
E. Write a frozen, GT-blind, costed two-axis comparison: proposal competition F3-native vs X-rational vs standard full-J; fixed-proposal evaluation competition exact-port vs ordinary goal-ROM vs A17 exchange vs full solver.
F. Give a CONDITIONAL GO / NO GO decision with exact missing evidence and smallest follow-up budget requested. Do not start paid physics merely because the archived 128/4096 maximum exists.

Output in a new non-overwriting research/a17x_f3_competition_audit/ directory:
- SOURCE_AND_TARGET_MANIFEST.md
- QUOTIENT_PORT_COMPATIBILITY.md
- COMPETITION_PROTOCOL_FROZEN.md
- ZERO_PHYSICS_ALGEBRA_CHECKS.json
- EVIDENCE_GAPS_AND_BUDGET.md
- GO_NO_GO.md
- START_HERE.md

If an authorized P1 follow-up begins later, start with a small pre-registered mechanism pilot. Allocate the same candidates to all evaluation arms; match material distance and total source/adjoint/setup cost, avoid leakage of full evaluators/GT into proposal selection, and insist on held-out independent material parents for any general claim. Posterior reference and NN require separate subsequent gates.
