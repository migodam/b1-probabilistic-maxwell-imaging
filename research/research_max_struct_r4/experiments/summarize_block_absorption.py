"""Read-only evidence summary; no physics, truth payloads or extra inversion."""
from pathlib import Path
from collections import Counter
import json, shutil, time

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]


def run():
    wall,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/block_absorption'
    load=lambda p:json.loads(p.read_text())
    block=load(out/'RECEIPT.json')
    aux=load(out/'AUXILIARY.json')
    check=load(out/'GRADIENT_WIRING.json')
    common=load(ROOT/'results/global_common_corrector/RECEIPT.json')
    r2=REPO/'research/research_max_struct_r2'
    r3=REPO/'research/research_max_struct_r3'
    case=block['configuration']['case_id']
    global_old=load(r2/f'results/blind/{case}/global/RECEIPT.json')
    global_aux=load(r2/f'results/blind/{case}/global/AUXILIARY.json')
    native=load(r2/f'results/blind/{case}/native/RECEIPT.json')
    scaled=load(r3/'results/development/scaled_original_restart/RECEIPT.json')
    rows=[]
    for label,rec,scope in [
        ('block_CIE_plus_common_corrector',block,'new full auxiliary plus corrector and wiring check'),
        ('R2_global_plus_common_corrector',common,'new corrector only; saved auxiliary reused'),
        ('published_R2_global',global_old,'historical full run'),
        ('published_R2_native',native,'historical full run'),
        ('R3_scaled_native_restart',scaled,'historical restart from published native endpoint')]:
        rows.append(dict(arm=label,original_Phi=rec['physical_objective'],
            physical_relative_residual=rec.get('physical_relative_residual',rec.get('physical_residual')),
            wall_seconds=rec['wall_seconds'],cpu_seconds=rec['cpu_seconds'],
            physical_LU=rec['physical_model_cost']['LU'],
            physical_complex_RHS=rec['physical_model_cost']['actual_complex_RHS'],cost_scope=scope))
    history=[x for x in aux['history'] if 'global_iteration' in x]
    counts=Counter(aux['execution_receipt']['operation_counts'])
    counts.update(check['operation_counts'])
    result=dict(status='PUBLIC_PHI_ONLY_NO_MATERIAL_SCORING',truth_access='NONE',rows=rows,
        stages=dict(block=aux['candidates'],R2_global=global_aux['candidates']),
        total_new_block_operation_counts_including_wiring=dict(counts),
        material_inner_stop_counts=dict(Counter(x['material_inner_message'] for x in history)),
        max_material_inner_iterations=max(x['material_inner_nit'] for x in history),
        maximum_current_step_objective_change=max(x['current_step_objective_change'] for x in history),
        maximum_material_step_objective_change=max(x['material_final_block_loss']-
            min(x['material_old_block_loss'],x['material_scalar_seed_block_loss']) for x in history),
        new_physical_cost_total=dict(LU=block['physical_model_cost']['LU']+common['physical_model_cost']['LU'],
            complex_RHS=block['physical_model_cost']['actual_complex_RHS']+common['physical_model_cost']['actual_complex_RHS']),
        qualification='single disclosed development case; unequal material-inner work; capped outer and final correctors; no material recovery or basin conclusion',
        additional_science_runs=0)
    (out/'COMPARISON.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    table='\n'.join(f"| {x['arm']} | {x['original_Phi']:.12g} | {x['physical_relative_residual']:.12g} | {x['physical_LU']} | {x['physical_complex_RHS']} |" for x in rows)
    note=f'''# Fixed block-loop absorption: implementation and recorded result

One configured 24×24 disclosed curvedbar development case was run, without reading truth. The new prediction and the global common-corrector control are sealed. No parameter tuning, CIE rerun for the control, formula identity suite, commit, or push was performed. Scientific conclusions and material scoring remain with the parent.

## Implementation and exact comparison boundary

`src/cie_block_absorption.py` retains the R2 global PR loop, fixed source metrics, native Fourier coefficient coordinates, DPCS offset, stage schedule, bounds, and original material prior. It imports the mature Fourier/DPCS helper and shared-real scalar CIE readout. Disjoint 2×2 spatial-cell groups give144 actual4×4 Green submatrices including self terms. These groups and their alignment are fixed.

The parent-supplied replacement is B=beta I-Gin, K=I+diag(m)B, R=beta K^-1(J-diag(m)Q). Current H/T use this same fixed-material transform. Material evaluation uses a=K^-H(R Omega), z=B R+beta Q, gradient=-rho Re sum(conj(a)*z)+prior gradient. The scalar CIE readout supplies a seed only; both previous and seed are scored with the exact block loss plus prior before L-BFGS-B(maxiter60), with fallback retained. This adopts known block absorption and CIE algebra and makes no originality assertion.

All three phases ran500 PR updates: mf3/beta1, mf6/current1/material0.2, mf10/current1/material0.2. DPCS rtol=1e-3, host=0.1, box=[0.02,20]. The prior uses the exact original `model.reg.evaluate` values/gradients; the unused dense Hessian built by the historical `model.regularizer` wrapper is omitted. Thus wall comparisons also include this implementation difference, and no speed claim is made.

The control loads the minimum-original-Phi R2 global auxiliary candidate (stage2, Phi0.011269442424001625), applies the identical R4 `cooperative_paths.refine(maximum=1000)`, and retains the published R2 global prediction as an incumbent. Both new paths use original Phi including prior for selection. Historical native and R3 scaled-native references have different cost histories and are context, not equal-compute arms.

## Public physical scores

| Arm | Original Phi including prior | Relative physical residual | Physical LU | Physical complex RHS |
|---|---:|---:|---:|---:|
{table}

Block auxiliary Phi was0.2645905236382217,0.030312610837575287,0.0007066885666717203 across stages; matched R2 global auxiliary Phi was0.2767745966664839,0.15300752495451386,0.011269442424001625. Both new physical correctors reached the1000-iteration cap. All three CIE outer stages reached500. Of1500 material inner solves,601 reached60 iterations,734 stopped on projected gradient, and165 on relative function reduction. No endpoint is asserted to be a certified local minimum or basin result.

## Validation, costs and provenance

The sole new gradient wiring check passed at beta0.2 with one fixed public-operator direction: relative error={check['relative_error']:.12g} at h=1e-6. Its three exact value/gradient calls cost864 small4×4 factorizations/solves and10368 complex RHS, with zero physical LU. There was no extra identity or error-floor suite.

Full block auxiliary work used61440 material value/gradient evaluations and53467 L-BFGS inner iterations. Including wiring: {counts['small4x4_batch_solve_calls']} batched calls, {counts['small4x4_factorizations']} small4×4 factorizations, {counts['small4x4_solves']} small4×4 solves, and {counts['small4x4_complex_RHS']} small-system complex RHS. NumPy's dense library solve refactors per call; this cost is explicitly separate from Maxwell LU/RHS. There were no inner fallbacks. Existing iteration receipts show decreasing fixed-current block material loss and decreasing fixed-material current objective at every recorded step; these are recorded run facts, not additional experiment tests.

Block process wall/CPU={block['wall_seconds']:.6f}/{block['cpu_seconds']:.6f}s; common-control new wall/CPU={common['wall_seconds']:.6f}/{common['cpu_seconds']:.6f}s. These include startup/imports and the artifacts produced inside each runner, with the summary postprocessing timed separately. Historical operator assembly was not rerun. New runs total2046 physical LU/48996 complex RHS. Both startup records imported successfully; failed_attempts=[] in both receipts. No costs for hypothetical unexecuted work are reported as measured.

`STARTUP.json`, `GRADIENT_WIRING.json`, `AUXILIARY_EXECUTION.json`, `AUXILIARY.json`, `CORRECTOR.json`, `RECEIPT.json`, `stage_0.npz` through `stage_2.npz`, `ALL_PHYSICAL_CANDIDATES.npz`, and `prediction.npz` preserve the evidence. `results/global_common_corrector/` holds its independent receipt and prediction. Executed source snapshots are saved; the current module differs from its executed copy only by explanatory docstring corrections. Frozen R2 CIE, mature baseline, and physical model source files were not modified.

Reproduce in a fresh output directory using the supplied `.venv_nn/bin/python`, CPU thread count1, and `experiments/block_absorption.py`, then `experiments/global_common_corrector.py`; both runners refuse to overwrite an existing disposition. `experiments/summarize_block_absorption.py` reads receipts only. No GT-bearing score files or truth payloads are read by any of these runners.
'''
    # Keep scientific variable notation intact while making prose readable.
    for a,b in [('give144','give 144'),('actual4×4','actual 4×4'),('ran500','ran 500'),
                ('stage2','stage 2'),('Phi0.','Phi 0.'),('was0.','was 0.'),
                ('the1000','the 1000'),('reached500','reached 500'),('Of1500','Of 1500'),
                ('solves,601','solves, 601'),('reached60','reached 60'),('iterations,734','iterations, 734'),
                ('and165','and 165'),('cost864','cost 864'),('and10368','and 10368'),
                ('used61440','used 61440'),('and53467','and 53467'),('total2046','total 2046'),
                ('count1','count 1')]:note=note.replace(a,b)
    (out/'METHOD_AND_RESULT.md').write_text(note)
    for name in ['block_absorption','global_common_corrector']:
        directory=ROOT/'results'/name
        shutil.copy2(ROOT/f'{name}_launch.log',directory/'RUN.log')
    snapshot=ROOT/'results/global_common_corrector/executed_source_snapshot'
    snapshot.mkdir(exist_ok=True)
    for p in [ROOT/'experiments/global_common_corrector.py',ROOT/'configs/global_common_corrector.json']:
        shutil.copy2(p,snapshot/p.name)
    post=dict(status='COMPLETED',additional_science_runs=0,truth_access='NONE',
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        physical_cost={},timer_scope='receipt summary and source/log copies; final postprocess receipt excluded')
    (out/'POSTPROCESS_RECEIPT.json').write_text(json.dumps(post,indent=2))
    print(json.dumps(dict(status=result['status'],rows=rows),indent=2))


if __name__=='__main__':run()
