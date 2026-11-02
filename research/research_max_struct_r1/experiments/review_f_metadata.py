"""Annotate saved receipts; no forward solve, unchanged scientific arrays."""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from tensor_families import rank_info
start=time.perf_counter();cpu=time.process_time()
O=ROOT/'results/families'
s=json.loads((O/'SUMMARY.json').read_text())
rankrows=[]
for r in s['rows']:
    n=r['mesh']['N'];d=np.load(O/(r['run_id']+'.npz'))
    checks={}
    for name,cols in [('one_old_source',[0]),('two_old_sources',[0,1]),('old_plus_independent',[0,1,2])]:
        checks[name]=rank_info(d['E'][:,cols],n)
    r['real_or_complex_rank']=checks['two_old_sources']
    rankrows.append(dict(run_id=r['run_id'],rank_cases=checks))
    r['E_state_dof_all_six_solved']=36*n
    r['same_field_annihilator_defect']=r.get('same_field_annihilator_defect',r['constitutive_defect']);r['constitutive_defect']=0.
    r['constitutive_defect_note']='j was constructed as chi E; original solver residual stored separately'
    leg=r['tensor_legal_margin']['branch']
    support=dict(real_positive=leg['real_min']>0,passive=leg['loss_min']>=0,nontrivial=r['branch_distance']['frobenius_rms']>0,
      algebraic_annihilator='exact by construction in real arithmetic; floating defects reported, not threshold-tuned')
    r['support_checks']=support
    r['status']='PASS' if all(support[k] for k in ['real_positive','passive','nontrivial']) else 'FAIL'
    r['cost']['wall_s_scope']='assembly+compute before array output; SUMMARY process total includes NPZ I/O'
    r['cost']['wall_including_npz_io_s']=None
    r['cost']['wall_including_npz_io_reason']='not separately timed in executed v1; process total includes I/O'
    r['validation']['H_relative_scope']='cell-center finite-DDA magnetic field, not continuum Hcurl'
    r['receipt_review']='v1.1 metadata-only correction; no numerical reexecution'
    (O/(r['run_id']+'.json')).write_text(json.dumps(r,indent=2))
r=s['interpolation'];r['same_field_annihilator_defect']=r.get('same_field_annihilator_defect',r['constitutive_defect']);r['constitutive_defect']=0.
r['transport']='coarse q interpolated/extrapolated, then analytic bump times q q^T; not direct tensor interpolation. Extrapolated centers outside bump support.'
r['receipt_review']='v1.1 metadata-only correction; no numerical reexecution'
(O/'F_interpolation.json').write_text(json.dumps(r,indent=2))
r=s['finite_noise'];r['unit_ray_excitation']=r.get('unit_ray_excitation',r['constitutive_defect']);r['constitutive_defect']=0.;r['state_dof']=6*r['mesh']['N']*3
q=r['noise_qualification'];leg=r['tensor_legal_margin']
r['support_checks']=dict(pair_bound_holds=q['pair_actual']<=q['pair_bound'],actual_observation_both_fit=q['both_fit'],strict_passivity=leg['loss_min']>0)
r['status']='PASS' if all(r['support_checks'].values()) else 'FAIL'
r['real_or_complex_rank']['smallest_returned_complex_singular_value']=r['real_or_complex_rank'].pop('smallest_nonzero_complex_singular_value',r['real_or_complex_rank'].get('smallest_returned_complex_singular_value'))
r['receipt_review']='v1.1 metadata-only correction; no numerical reexecution'
(O/'F_noise_ray.json').write_text(json.dumps(r,indent=2))
(O/'SUMMARY.json').write_text(json.dumps(s,indent=2))
(O/'RANK_CLASSES.json').write_text(json.dumps(rankrows,indent=2))
receipt=dict(operation='saved-array rank class analysis and receipt metadata review',new_forward_LU=0,new_forward_RHS=0,local_small_SVD=6*sum(r['mesh']['N'] for r in s['rows'])+24,
             wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,stage='scoring',
             changes=['correct E DOF','separate constitutive residual and VU excitation','legality/bound support booleans','rank label/zero-row edge case','interpolation and timer scope'])
(ROOT/'audit/F_METADATA_REVIEW.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt))
