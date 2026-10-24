# Paired actual-physics refinement: inert preregistration

Status: **READY BUT UNRELEASED**. `configs/PHYSICS_TRANSPORT_PAIRED_V1.json` has `released:false`. No fit, optimizer update, child process or physical query was launched during preparation. The parent owns release and every scientific interpretation.

The serial controller is `tools/queue_r4_physics_transport.py`. It uses the existing trainer unchanged, with its current exact full66 objective, native original complex128 observations, immediate current/field/query persistence, original TRAIN normalizers, inclusive NLL batch128, physical microbatch4,3000 updates, LR5e-5, gradient clip5, fixed lambda schedule, validation at0/300/.../3000 and fixed NLL+.1 reverse-objective selection. Every fit has a2400s Book quota; physics/control RHS hard caps are256000/20000. Planned complete costs are203264/11264 per pair arm. `stored` codec is frozen for both arms. The controller does not wrap a trainer in an analysis Book; each child owns its standalone physics Book and both physical/training locks.

## Frozen families and order

| Family | TRAIN / VAL | Parameters | Original checkpoints |
|---|---:|---:|---|
| raw spline |32768 /2048 |682860 |`density_spline_32768_s202610250{1,2,3}` |
| matched raw diagonal Gaussian |32768 /2048 |682943 |`extended_diag_gaussian_s202610250{1,2,3}` |
| anchor spline, coordinate none |1024 /256 |691308 |`anchor_spline_none_s202610250{1,2,3}` |

All nine original source fits are CLOSED COMPLETE; registry records checkpoint/result/receipt byte+mtime identities and exact original seed. No source is chosen by new outcomes. The raw-y and anchor data panels/information sets remain distinct; this plan is not an equal-data comparison between them.

First run seed2026102501 pairs in the table's order, physics then supervised within each pair. **Finish all three seed1 pairs before any replication.** Compute each gate only from the physics fit's existing trainer-selected checkpoint and its initialization; retain its paired supervised fit as the comparison. Do not choose another checkpoint to clear the gate. Both seed1 arms must be COMPLETE with3000 updates.

A family is eligible iff all three conditions hold:

1. `(initial_composite-selected_composite)/max(abs(initial_composite),1e-12) >= .05`.
2. `selected_meanPhi <= .8 * initial_meanPhi`.
3. `selected_proper_NLL <= initial_proper_NLL +10` nats in raw Gaussianized dx66.

For eligible families, run both original-source seed2026102502 and2026102503 pairs in family order, physics then supervised. Seed2 outcome does not determine whether seed3 runs. Ineligible families get no repeated fits. Gates use the already exposed fixed validation panel; they are continuation decisions, not fresh holdout or posterior-qualification evidence.

Ceiling:9 physics and9 supervised fits total, at most1930752 planned RHS,2484000 aggregate hard-cap RHS, and43200 aggregate per-fit wall quota seconds. These are ceilings, not instructions to spend unused resources. A failed/open/inconsistent fit stops the controller with HOLD; there is no automatic replay, replacement, hyperparameter revision or failed-slot deletion.

## Release and execution

The registry freezes the unchanged trainer and supporting source byte/mtime identities, target config, same-input VJP audit, original1056 observation identities, model inputs/outputs, common protocol and gate before any fit. The original64-RHS real preflight FAIL is preserved. The separate same-input correction audit is `PASS_CORRECTED_REFERENCE`, using saved paid VJP with zero extra physics and no tolerance relaxation. This gate permits implementation use; it does not qualify a posterior.

Review the registry, then the parent must explicitly set `released:true`. Default invocation only validates metadata and prints readiness:

```text
<research-python> tools/queue_r4_physics_transport.py
```

Only the explicit execution command launches the released serial plan:

```text
<research-python> tools/queue_r4_physics_transport.py --execute
```

Do not place it under `run_r4_job.py`. It waits for both campaign slots without holding them; the child reacquires its required locks and a race/failure stops the controller. Logs, per-fit frozen releases, gate records, progress and terminal controller result go under `results/physics_transport_paired_v1`; outputs go under `models/physics_transport/paired_v1`. Idle controller wall is recorded separately and not summed again as child fit wall. Source/registry changes after start fail closed. An explicit `--resume` may skip only verified COMPLETE jobs with matching frozen releases; failed/open jobs never restart automatically.

## Readiness evidence

`AUDIT/physics_transport_campaign_readiness_1792877680066764000/READINESS.json`: PASS_READY_INERT. Analysis-only Book I COMPLETE, wall1.854s, CPU1.735s, peakRSS227606528B. Eighteen conditional plans match all caps. Mock gate tests cover all-pass, each separate failure, exact thresholds, failed-fit rejection and rejection of alternate-checkpoint reselection. Unreleased execution was rejected before any child launch. Actual Target instances0, physical RHS0, optimizer updates0. The trainer/target/math and old reports were not changed.
