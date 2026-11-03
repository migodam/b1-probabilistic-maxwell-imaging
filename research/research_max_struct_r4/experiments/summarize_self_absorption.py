"""Summarize the fixed self-only control without truth reads or physics calls."""
from pathlib import Path
from collections import Counter
import ast,json,shutil,time

ROOT=Path(__file__).resolve().parents[1]


def run():
    wall,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/self_absorption'
    load=lambda p:json.loads(p.read_text())
    receipt=load(out/'RECEIPT.json');aux=load(out/'AUXILIARY.json')
    block=load(ROOT/'results/block_absorption/RECEIPT.json')
    control=load(ROOT/'results/global_common_corrector/RECEIPT.json')
    counts=aux['execution_receipt']['operation_counts']
    history=[x for x in aux['history']if 'global_iteration'in x]
    old_source=(ROOT/'results/block_absorption/executed_source_snapshot/cie_block_absorption.py').read_text()
    new_source=(out/'executed_source_snapshot/cie_block_absorption.py').read_text()
    def function(source,name):
        tree=ast.parse(source)
        node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name==name)
        return ast.dump(node,include_attributes=False)
    shared=['block_material_value_gradient','block_material_readout']
    unchanged={name:function(old_source,name)==function(new_source,name)for name in shared}
    fields=['stages','host','material_box','dpcs_rtol','coefficient_tol','miniter',
            'material_inner_maxiter','material_inner_ftol','material_inner_gtol','prior',
            'source_metric_policy','corrector','corrector_maxiter']
    configuration_matches={key:receipt['configuration'][key]==block['configuration'][key] for key in fields}
    rows=[dict(arm=name,original_Phi=r['physical_objective'],
        relative_physical_residual=r['physical_relative_residual'])
        for name,r in [('self_only',receipt),('block4',block),('R2_global_common_corrector',control)]]
    result=dict(status='SEALED_PUBLIC_PHYSICAL_COMPARISON_NO_GT_SCORE',truth_access='NONE',
        rows=rows,configuration_matches=configuration_matches,
        shared_material_functions_AST_unchanged=unchanged,
        material_inner_stops=dict(Counter(x['material_inner_message']for x in history)),
        small4x4_factorizations=counts.get('small4x4_factorizations',0),
        small4x4_solves=counts.get('small4x4_solves',0),
        scalar_complex_divisions=counts['scalar_complex_divisions'],
        note='This one self-only control did not reproduce block4 physical fit; no general necessity, uniqueness, basin or material-recovery conclusion.',
        additional_numerical_checks=0,additional_science_runs=0)
    (out/'COMPARISON.json').write_text(json.dumps(result,indent=2))
    table='\n'.join(f"| {r['arm']} | {r['original_Phi']:.12g} | {r['relative_physical_residual']:.12g} |" for r in rows)
    note=f'''# Self-only absorption control

The single frozen self-only run is complete and its prediction is sealed. No truth payload, material score, parameter tuning, additional derivative/identity suite, commit, or push was used.

The scientific change is `gin_mode='self'`: Gin contains only diag(G), with all intercell entries removed. The operator applies P_i=beta/[1+(beta-G_ii)m_i] by direct complex division. Its adjoint is also diagonal division. These are charged as actual scalar divisions, with zero small4×4 factorizations or solves. Default `gin_mode='block4'` retains the previous block branch. The original block execution source, frozen configuration, receipts and prediction remain untouched.

The generic exact material value/gradient and mature scalar CIE seed plus L-BFGS readout functions are AST-identical to the original executed block implementation. The proposed exact rational self readout was not implemented. All matched scientific fields listed in `COMPARISON.json` are unchanged: source metrics, prior, bounds, initialization, Fourier/DPCS, three500-update stages, max60 material iterations, and common scaled original-Phi corrector(maxiter1000). Candidate selection uses original Phi including prior and retains the original host incumbent. The one prior block-gradient wiring receipt is referenced; no new numerical check was run.

| Arm | Original Phi including prior | Relative physical residual |
|---|---:|---:|
{table}

Self auxiliary stage Phi values were0.26587463551255636,0.06756084411575954,0.010801894717538736. All three stages reached500 updates. The common corrector reached1000 iterations with1048 evaluations. Among1500 material solves,676 hit the60-iteration cap,711 stopped on projected gradient and113 on relative function reduction; zero fallbacks occurred. These caps are preserved, not convergence claims.

This run did not reproduce the block4 physical-fit improvement using only the self term. It does not establish that intercell blocks are universally necessary or prove a basin or material-recovery result. Parent scientific scoring and conclusions are separate.

Measured process wall/CPU: {receipt['wall_seconds']:.6f}/{receipt['cpu_seconds']:.6f}s. Physical Maxwell cost: {receipt['physical_model_cost']['LU']} LU and {receipt['physical_model_cost']['actual_complex_RHS']} complex RHS. Self diagonal work: {counts['self_diagonal_solve_calls']} application calls and {counts['scalar_complex_divisions']} actual complex scalar divisions ({counts['self_forward_scalar_divisions']} forward, {counts['self_adjoint_scalar_divisions']} adjoint). Material evaluations: {counts['material_inner_evaluations']}; inner L-BFGS iterations: {counts['material_inner_LBFGS_iterations']}. Source-metric Cholesky and receiver/APCS SVD counts remain separately reported. Startup imported successfully and failed_attempts=[]; operator assembly was not rerun.

Files include `FROZEN_CONFIG.json`, `STARTUP.json`, `AUXILIARY_EXECUTION.json`, `AUXILIARY.json`, `CORRECTOR.json`, `RECEIPT.json`, three stage-current NPZs, `SELECTED_BEFORE_CORRECTOR.npz`, `ALL_PHYSICAL_CANDIDATES.npz`, `prediction.npz`, `RUN.log` and `executed_source_snapshot/`. The original block note predates this extension; its executed snapshot is authoritative for the old run, while this run has its own complete source snapshot.

Command: the scientific `.venv_nn/bin/python` with CPU thread count1 runs `experiments/block_absorption.py --variant self_absorption`. An existing result receipt prevents overwrite. `experiments/summarize_self_absorption.py` reads saved receipts and source only, and its postprocessing timing is separate.
'''
    for a,b in [('three500','three 500'),('max60','max 60'),('were0.','were 0.'),
                ('reached500','reached 500'),('reached1000','reached 1000'),('with1048','with 1048'),
                ('Among1500','Among 1500'),('solves,676','solves, 676'),('the60','the 60'),
                ('cap,711','cap, 711'),('and113','and 113'),('count1','count 1')]:note=note.replace(a,b)
    (out/'METHOD_AND_RESULT.md').write_text(note)
    shutil.copy2(ROOT/'self_absorption_launch.log',out/'RUN.log')
    post=dict(status='COMPLETED',truth_access='NONE',additional_science_runs=0,
        physical_cost={},wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
    (out/'POSTPROCESS_RECEIPT.json').write_text(json.dumps(post,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
