"""Do not overwrite the original failed F64 intervals; append bounded audit.

Uses M4's own sealed query/potential transcript, never truth or a physics
teacher. All certificate computation must count toward deployment if used.
"""
from pathlib import Path
from decimal import Decimal
import argparse,json,time
import numpy as np
from bgga.decimal_bounds import normalized_intervals,independent_reference
from g1_replay import audited_inputs, _load_keys
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--old-bank-root',type=Path,default=None,
        help='explicit relocated old finite_bank/online directory; no automatic data search')
    args=parser.parse_args()
    source=ROOT/'results/g1_replay_001'
    assert json.loads((source/'ALL_ONLINE_OUTPUTS_SEALED.json').read_text())['all_method_outputs_written']
    out=ROOT/'results/g1_normalization_repair';out.mkdir(parents=True,exist_ok=False)
    inputs=audited_inputs(json.loads((ROOT/'research/delegated/source_audit/CACHE_MAP.json').read_text()),args.old_bank_root)
    grid=inputs['FINITE_PRIOR_GRID']
    p=_load_keys(grid['path'],('prior_weight',))['prior_weight']
    rows=[]
    for parent in range(12):
        for rep in range(2):
            t=time.process_time();case=source/f'parent_{parent:02d}_repeat_{rep:02d}'
            with np.load(case/'M4.npz',allow_pickle=False) as z:
                order=z['query_indices'];query_ll=z['query_log_likelihood'];qpoint=z['joint_weights']
            ll=np.full(len(p),np.nan);ll[order]=query_ll
            if np.any(~np.isfinite(ll)):raise ValueError('This closeout audit expects the complete original trace, no added physical queries')
            reference=independent_reference(p,ll)
            max_point_error=max(abs(Decimal.from_float(float(q))-r) for q,r in zip(qpoint,reference))
            checkpoint_rows=[];arrays={}
            for n in [16,64,256,1024,2048,2949,2950]:
                known=np.zeros(len(p),bool);known[order[:n]]=True
                bounds=normalized_intervals(p,ll,known)
                lower_bad=[];upper_bad=[]
                for i,(lo,hi,r) in enumerate(zip(bounds['lower'],bounds['upper'],reference)):
                    if Decimal.from_float(float(lo))>r:lower_bad.append(i)
                    if r>Decimal.from_float(float(hi)):upper_bad.append(i)
                checkpoint_rows.append(dict(query_count=n,lower_violations=len(lower_bad),upper_violations=len(upper_bad),
                    unknown_mass_upper=bounds['unknown_mass_upper'],scope=bounds['label']))
                arrays[f'lower_{n}']=bounds['lower'];arrays[f'upper_{n}']=bounds['upper']
            np.savez_compressed(out/f'parent_{parent:02d}_repeat_{rep:02d}.npz',**arrays)
            rows.append(dict(parent=parent,repeat=rep,checkpoints=checkpoint_rows,
                maximum_original_point_error_to_80_digit_normalization=float(max_point_error),
                original_queries=int(len(order)),new_physics_queries=0,
                measured_certificate_and_independent_audit_cpu_seconds=time.process_time()-t))
            (out/'RESULTS.json').write_text(json.dumps(rows,indent=2)+'\n')
            print(json.dumps(rows[-1]),flush=True)
    result=dict(status='NORMALIZATION_AUDIT_COMPLETE',rows=24,
        any_corrected_bound_violation=any(c['lower_violations'] or c['upper_violations'] for r in rows for c in r['checkpoints']),
        old_raw_violations_retained=True,old_results_modified=False,new_FE=0,new_RHS=0,
        scope='conditional finite numerical log-potential model only',
        not_certified=['Maxwell/FEM forward error','formation of float64 Gaussian potential','continuum PDE'],
        deployment_attribution='If used online, the normalization cost is additional M4 certificate cost, never free offline preparation',
        cost_route_still_no_go='All original2950 queries per case still required, including after repair')
    (out/'DECISION.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
