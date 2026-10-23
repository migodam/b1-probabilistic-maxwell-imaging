# W execution integrity audit

Status: **FINAL_REPLAY_WITH_DECLARED_LIMITATIONS**. Snapshot: 2026-10-23T18:43:18.475043+00:00.

Accounting totals below come only from the root ledger. Receipt CPU is not charged a second time. This report does not judge scientific acceptance.

| Stage | RHS / cap | LU / cap | settled all-process CPU / cap (s) | sealed |
|---|---:|---:|---:|---|
| R | 12000 / 12000 | 1195 / 3000 | 2251.516 / 14400 | True |
| B2 | 16140 / 24000 | 1695 / 6000 | 8386.091 / 28800 | False |
| C | 2700 / 8000 | 360 / 2000 | 1992.751 / 10800 | False |
| P1A | 13880 / 24000 | 3470 / 6000 | 15283.426 / 28800 | False |
| P1B | 15640 / 24000 | 3491 / 6000 | 17767.982 / 28800 | False |

| B2 scope | RHS / cap | LU (no separate scope cap declared) |
|---|---:|---:|
| cal_2026102202_narrow | 5220 / 6000 | 510 |
| joint_2026102202_narrow | 2420 / 3500 | 155 |
| cavity_2026102202_narrow | 2420 / 3500 | 155 |
| cal_2026102201_narrow | 1240 / 5000 | 115 |
| joint_2026102201_narrow | 2440 / 3500 | 160 |
| validation | 2400 / 2500 | 600 |

## Recorded checks

- Reservation ID duplicates: 0.
- Solver events matched to prior reservations: 22431.
- Solver events without prior reservations: 0.
- Maximum logged RHS-relative solver residual: 2.06632e-12 at ledger line 71552.
- Indexed endpoint states: {'DECLARED_COMPLETE_FILES_PRESENT': 1893}.
- Safe metadata sample: 24; arrays opened: 0.
- Root settled CPU: 45681.764441s; nested solver CPU: INFORMATIONAL_NONADDITIVE.

## Findings and pending settlement

- **UNVERIFIED / STAGES_NOT_SEALED**: Stage totals are a live/unsealed snapshot, not final audited expenditure

## Failed runs remain paid

- `p1b-independent-offline-risk`: FAILED, CPU 0.06204s, contextual RHS 0 / LU 0, MATCHED_ROOT_LEDGER; ledger lines [70135].
- `r0-truth-closure`: FAILED_PRE_EXEC, CPU 0.043587s, contextual RHS 0 / LU 0, MATCHED_ROOT_LEDGER; ledger lines [1].
- `r3-graded-optimizer-controls`: FAILED, CPU 69.72585099999999s, contextual RHS 772 / LU 74, MATCHED_ROOT_LEDGER; ledger lines [8708].

## Concrete current-code evidence

- reserves_are_append_only: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1/src/wr/accounting.py` (CURRENT_SOURCE_ONLY).
  - Line 36: `def reserve(self,phase,label,rhs=0,factorizations=0,**meta):`
  - Line 48: `row=dict(event='budget_reject' if reason else 'reserve',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),stage=self.stage,phase=phase,label=label,rhs=int(rhs),factorizations=int(factorizations),id=uuid.uuid4().hex,**meta)`
  - Line 50: `f.seek(0,2); f.write(json.dumps(row,default=encode,allow_nan=False)+'\n'); f.flush()`
- ledger_is_authoritative: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1/src/wr/accounting.py` (CURRENT_SOURCE_ONLY).
  - Line 23: `def totals(self):`
  - Line 29: `result[stage]={k:sum(r.get(k,0) for r in rr if r.get('event')=='reserve') for k in ['rhs','factorizations']}`
  - Line 30: `result[stage]['cpu_seconds']=sum(r.get('charged_cpu_seconds',0) for r in rr)`
- scope_cap_is_not_reset: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1/src/wr/accounting.py` (CURRENT_SOURCE_ONLY).
  - Line 54: `class ScopedBook(Book):`
  - Line 65: `raise BudgetStop('scope_rhs_limit: '+self.scope)`
- all_child_process_CPU_and_single_ledger_charge: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1/tools/run_job.py` (CURRENT_SOURCE_ONLY).
  - Line 17: `start=time.perf_counter(); before=resource.getrusage(resource.RUSAGE_CHILDREN)`
  - Line 41: `cpu=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime`
  - Line 43: `write(dest/'RECEIPT.json',row); row.pop('stage'); book.event(**row)`
- frequency_files_before_endpoint_completion: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1/src/wr/physics.py` (CURRENT_SOURCE_ONLY).
  - Line 100: `np.savez_compressed(dest/f'f{idx}.npz',x=x,frequency=f,data=s.data,probe_data=s.probe_data,epsilon=s.epsilon,points=s.material_points,weights=s.material_weights,field=s.field,current=s.current,scattered_dofs=s.scattered_dofs,mesh_p=p.mesh.p,mesh_t=p.mesh.t,residuals=s.residuals)`
  - Line 111: `write(dest/'COMPLETE.json',dict(identity=self.identity,frequency_count=len(self.freq),jacobian=jac,totals=self.book.totals()))`
- reserve_before_action_and_nested_wall_timers: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_continuation/src/wq/fem.py` (CURRENT_SOURCE_ONLY).
  - Line 550: `self.book.reserve(self.phase, 'fem_forward_sources', rhs=count,`
  - Line 479: `factor = splu(matrix, permc_spec='COLAMD')`
  - Line 591: `timings = dict(setup_seconds=self.setup_seconds,`
  - Line 602: `self._event('backend_solve', freq=self.freq, rhs=count,`

## Coverage limits

- Root ledger reserve rows, including failed attempts, are authoritative. Rejections are not successful reservations.
- Current parse/ID uniqueness and append-mode source evidence do not prove historical append-only immutability; no hash or trusted-prefix comparison performed.
- Wrapper receipts have no child PID/reserve-ID foreign key; child attribution is stage/time-span containment and may remain ambiguous.
- Endpoint existence/nonempty checks are not array decompression, numerical content validation or scientific correctness.
- Only ledger-indexed Forward endpoints are inventoried. Legacy unindexed endpoint persistence is NOT exhaustively verified.
- Private generated-data endpoint metadata, OFFLINE_TRUTH, GENERATION_PLAN, observations and job logs are never opened.
- Live/unreceipted process CPU is not yet settled in root ledger; no assertion that charged remaining CPU equals live headroom.
- Nested backend timers overlap, usually measure wall time, and must not be added to all-process CPU charges.
- Current source snippets document semantics; historical source snapshots are checked for existence, not equated to current code.

Audit actions: zero RHS, zero LU, zero GT files, no hashes; no physics/configuration/ledger/job/array changes.

Replay with `python tools/audit_w_execution.py --final` after settlement. Live stages are not failed merely because they are unfinished.
