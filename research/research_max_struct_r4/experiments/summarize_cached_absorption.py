"""Read saved execution evidence only; no truth reads or new physics calls."""
from pathlib import Path
from collections import Counter
import ast,json,shutil,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def run():
    wall,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/block_absorption_cached';old=ROOT/'results/block_absorption'
    load=lambda p:json.loads(p.read_text())
    new_rec,old_rec=load(out/'RECEIPT.json'),load(old/'RECEIPT.json')
    new_aux,old_aux=load(out/'AUXILIARY.json'),load(old/'AUXILIARY.json')
    new_corr,old_corr=load(out/'CORRECTOR.json'),load(old/'CORRECTOR.json')
    new_counts=new_aux['execution_receipt']['operation_counts']
    old_counts=old_aux['execution_receipt']['operation_counts']
    rows=[]
    for name,rec,aux,corr in [('direct',old_rec,old_aux,old_corr),
                              ('cached_inverse',new_rec,new_aux,new_corr)]:
        count=aux['execution_receipt']['operation_counts']
        rows.append(dict(arm=name,original_Phi=rec['physical_objective'],
            relative_physical_residual=rec['physical_relative_residual'],
            full_linear_relative_residual=rec['linear_relative_residual'],
            process_wall_seconds=rec['wall_seconds'],process_cpu_seconds=rec['cpu_seconds'],
            auxiliary_wall_seconds=aux['execution_receipt']['wall_seconds'],
            auxiliary_cpu_seconds=aux['execution_receipt']['cpu_seconds'],
            corrector_wall_seconds=corr['wall_seconds'],corrector_evaluations=corr['evaluations'],
            physical_LU=rec['physical_model_cost']['LU'],
            physical_complex_RHS=rec['physical_model_cost']['actual_complex_RHS'],
            auxiliary_small4x4_factorizations=count['small4x4_factorizations'],
            auxiliary_small4x4_solve_RHS=count['small4x4_complex_RHS'],
            auxiliary_inverse_matrix_apply_RHS=count.get('small4x4_matrix_apply_complex_RHS',0),
            material_evaluations=count['material_inner_evaluations'],
            material_LBFGS_iterations=count['material_inner_LBFGS_iterations']))
    with np.load(old/'prediction.npz',allow_pickle=False)as f:old_m=f['m'].copy()
    with np.load(out/'prediction.npz',allow_pickle=False)as f:new_m=f['m'].copy()
    fields=['stages','host','material_box','dpcs_rtol','coefficient_tol','miniter',
            'material_inner_maxiter','material_inner_ftol','material_inner_gtol','prior',
            'source_metric_policy','corrector','corrector_maxiter','material_seed']
    matches={k:old_rec['configuration'][k]==new_rec['configuration'][k]for k in fields}
    def derivative(source):
        node=next(x for x in ast.parse(source).body if isinstance(x,ast.FunctionDef)
                  and x.name=='block_material_value_gradient')
        return ast.dump(node,include_attributes=False)
    derivative_same=derivative((old/'executed_source_snapshot/cie_block_absorption.py').read_text())==derivative(
        (out/'executed_source_snapshot/cie_block_absorption.py').read_text())
    history=[x for x in new_aux['history']if 'global_iteration'in x]
    result=dict(status='SEALED_IMPLEMENTATION_REPLAY_NO_GT_SCORE',truth_access='NONE',rows=rows,
        scientific_configuration_matches=matches,material_derivative_AST_unchanged=derivative_same,
        full_run_observed_wall_reduction_fraction=1-new_rec['wall_seconds']/old_rec['wall_seconds'],
        auxiliary_observed_wall_reduction_fraction=1-new_aux['execution_receipt']['wall_seconds']/old_aux['execution_receipt']['wall_seconds'],
        material_difference_from_direct=dict(l2=float(np.linalg.norm(new_m-old_m)),
            relative_l2=float(np.linalg.norm(new_m-old_m)/np.linalg.norm(old_m)),
            max_abs=float(np.max(np.abs(new_m-old_m))),semantics='prediction difference, not ground-truth error'),
        cache_statistics={k:v for k,v in new_counts.items()if 'cache' in k or 'resident' in k or 'replacements'in k},
        material_inner_stops=dict(Counter(x['material_inner_message']for x in history)),
        efficiency_assessment='Observed lower auxiliary time and fewer actual factorizations/RHS in one replay; no universal speed claim.',
        attribution_limit='Unchanged physical corrector wall also changed; trajectory and evaluation counts differ; no repeated controlled benchmark.',
        old_wiring_check_wall_seconds=old_rec['gradient_check']['wall_seconds'],
        new_derivative_suite='NOT_RUN; unchanged derivative',additional_science_runs=0,
        original_direct_and_self_outputs_preserved=True)
    (out/'COMPARISON.json').write_text(json.dumps(result,indent=2))
    table='\n'.join(f"| {x['arm']} | {x['original_Phi']:.12g} | {x['relative_physical_residual']:.12g} | {x['auxiliary_wall_seconds']:.6f} | {x['process_wall_seconds']:.6f} |"for x in rows)
    note=f'''# Cached inverse implementation replay

Exactly one authorized development replay was completed. No truth payloads or material scores were read; no numerical derivative suite, parameter tuning, repeated timing experiment, publication or commit was performed. Previous direct-block and self-only results and executed-source snapshots are preserved.

## Implementation scope

`block_solver='cached_inverse'` is explicit; the default remains `direct`. For a cache miss, the code calls the existing NumPy batched `solve(K,I4)` on each actual4×4 block. It retains B,K,inverse and a copy of m in an LRU containing at most2 beta entries, one latest exact material per beta. Hits require equal beta and `np.array_equal(m)`. Residual/adjoint and H/T applications use inverse matrix products, with the adjoint applying inverse^H. No new factorization algorithm, LU kernel, precision, spatial block, or physics backend was introduced.

The optional value/gradient memo stores at most2 exact material points within each fixed-J/Q readout and resets before the next readout. It only reuses exact repeated calls; gradients are copied when returned from cache. The mathematical derivative function is AST-identical to the original executed source. Seed, prior, bounds, source metrics, tolerances,1500 outer updates and max1000 original-Phi corrector remain unchanged. No explicit rounding or quantization is added. Inverse multiplication changes floating-point operation order, so the trajectory and endpoints are retained separately.

## Measured physical outcome and timing

| Solver | Original Phi including prior | Relative physical residual | Auxiliary wall seconds | Process wall seconds |
|---|---:|---:|---:|---:|
{table}

Observed auxiliary wall reduction was {100*result['auxiliary_observed_wall_reduction_fraction']:.3f}%; process wall reduction was {100*result['full_run_observed_wall_reduction_fraction']:.3f}%. This is one replay, not a general speed guarantee. In particular, the unchanged physical corrector itself changed from{old_corr['wall_seconds']:.6f}s/{old_corr['evaluations']} evaluations to{new_corr['wall_seconds']:.6f}s/{new_corr['evaluations']} evaluations. Consequently the entire time reduction cannot be attributed to cached auxiliary kernels. The old process also included its{result['old_wiring_check_wall_seconds']:.6f}s wiring check, which was not repeated.

The cached final full linear-system relative residual is{new_rec['linear_relative_residual']:.12g}. The final prediction differs from direct by relative material L2={result['material_difference_from_direct']['relative_l2']:.12g} and max absolute material difference={result['material_difference_from_direct']['max_abs']:.12g}; these compare two predictions, not truth. Both final correctors reached1000 iterations, all three outer stages reached500, and the cached run had597 inner solves reach60 iterations. The physical-fit scale is retained, but no GT-recovery equivalence or superiority is asserted.

## Actual costs

Direct auxiliary work used{old_counts['small4x4_factorizations']} small4×4 factorizations and{old_counts['small4x4_complex_RHS']} factor-solve RHS. Cached work used{new_counts['small4x4_factorizations']} factorizations/solves with exactly4 RHS per block, totaling{new_counts['small4x4_complex_RHS']} inverse-solve RHS. Its additional{new_counts['small4x4_matrix_apply_complex_RHS']} matrix-apply RHS are separately counted and are not triangular solves. The inverse cache had{new_counts['small4x4_inverse_cache_hits']} exact hits,{new_counts['small4x4_inverse_cache_misses']} misses, and max resident entries={new_counts['small4x4_inverse_max_resident_beta_entries']}. Value/gradient memoization served{new_counts['material_value_gradient_cache_hits']} exact repeat hits; actual material evaluations={new_counts['material_inner_evaluations']}, inner L-BFGS iterations={new_counts['material_inner_LBFGS_iterations']}. No inner fallback occurred.

The new full run charged{new_rec['physical_model_cost']['LU']} physical Maxwell LU and{new_rec['physical_model_cost']['actual_complex_RHS']} physical complex RHS. Process CPU={new_rec['cpu_seconds']:.6f}s; auxiliary CPU={new_aux['execution_receipt']['cpu_seconds']:.6f}s. Startup succeeded; failed_attempts=[]; frozen operator assembly was not repeated. Snapshot and summary postprocessing have separate nonphysics receipts. Prior original costs remain recorded, rather than overwritten by the faster observation.

`RECEIPT.json`, `AUXILIARY_EXECUTION.json`, `AUXILIARY.json`, `CORRECTOR.json`, `FROZEN_CONFIG.json`, `STARTUP.json`, `RUN.log`, three stage NPZs, candidate NPZs, `prediction.npz`, and `executed_source_snapshot/` preserve the run. `COMPARISON.json` includes the configuration/source isolation evidence and all measured counts. Use the scientific Python with CPU threads1 to run `experiments/block_absorption.py --variant block_absorption_cached` in a fresh output directory; the runner refuses overwrite.
'''
    for a,b in [('actual4×4','actual 4×4'),('most2','most 2'),('tolerances,1500','tolerances, 1500'),
        ('max1000','max 1000'),('from'+str(''), 'from')]:
        note=note.replace(a,b)
    # Add spaces around recorded counts without altering numeric content.
    import re
    note=re.sub(r'(from|to|its|is|reached|had|reach|used|and|totaling|additional|hits,|served|charged|threads)(?=[0-9])',r'\1 ',note)
    (out/'METHOD_AND_RESULT.md').write_text(note)
    shutil.copy2(ROOT/'block_absorption_cached_launch.log',out/'RUN.log')
    post=dict(status='COMPLETED',truth_access='NONE',additional_science_runs=0,physical_cost={},
        wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
    (out/'POSTPROCESS_RECEIPT.json').write_text(json.dumps(post,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
