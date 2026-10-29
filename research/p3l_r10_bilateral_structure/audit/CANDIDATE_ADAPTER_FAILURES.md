# Candidate adapter failure disclosure

Two integration failures occurred in the fresh adapter, and both were repaired without changing methods, seeds, statistics, or scientific stopping policies.

1. The first generator process failed in its first bookkeeping event because `phase` was provided twice. The ledger was empty, no truth state existed, and no forward, adjoint, or LU was reserved. The empty attempt and source/config revisions are preserved under `candidate_fresh_prephysics_bookkeeping_fix/`.
2. The initial online fresh process omitted `ph.y=y`, which the historical R9 corrector reads internally. It stopped after the parent000 shared anchor, all three new-arm K0 endpoints, and the first old-arm raw/J pair. Cumulative cost at stop was2288 RHS256 and86 LU256, including data generation. The exact arrays, ledger, original log, and source/config revision are preserved. Resume phase `ONLINE_RESUME1` reuses the completed anchor and three completed new proposal records, sets `ph.y`, and replays the interrupted old candidate from its frozen initial point. The replay charges an additional40 RHS256 and1 LU256. Those wasted original raw/J costs remain in the global ledger; they are additional adapter overhead rather than hidden inside the replay's six-forward/four-J method allowance. No cached raw-state injection was used.

The phase label, original role fields, cumulative wall/CPU accounting, required historical backend `y` field, counter keys, and distinct endpoint names were inspected in code. No microscopic physics identity test was added. No fresh private truth was read during these repairs. The online prediction seal is required before the separate scorer can read private truth.

`failed_actions=0` in the physics resource counter means no failed physical solve; it does not mean that these two software integration failures did not occur. The two failures are part of the experimental audit history and must remain visible in delivery.
