"""Prepare compact readable endpoints; preserve raw runs in the full release ZIP."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    manifests=sorted((ROOT/'results').glob('*/MANIFEST.json'))
    count=0
    for path in manifests:
        manifest=json.loads(path.read_text())
        for row in manifest['rows']:
            source=ROOT/row['output'];result=json.loads(source.read_text())
            for field in ('trajectory','material_trajectory'):result.pop(field,None)
            result['raw_archive_path']=str(source.relative_to(ROOT))
            result['full_material_trajectory_archive_path']=str(source.with_suffix('.trajectory.npz').relative_to(ROOT))
            result['publication_note']='Raw trials and full trajectories retained in FULL_REPRODUCIBILITY.zip; compact endpoints contain no recomputed physics.'
            dest=ROOT/'published_results'/path.parent.name/source.name
            dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(result,indent=2)+'\n')
            count+=1
    runs=[]
    for path in sorted((ROOT/'results').glob('*/SCORES.json')):
        run=json.loads(path.read_text());run['run_name']=path.parent.name;runs.append(run)
    fresh=next(r for r in runs if r['run_name']=='fresh_v3')
    summary=dict(schema='r4-final-experiment-results',
       status='HOLD_FOR_PAPER; NO_FRESH_GT_IMPROVEMENT',
       scope='29 exposed-truth development runs and 18 frozen fresh predictions, 576 real materials, unchanged nonlinear VIE',
       runs=runs,rows=[row for run in runs for row in run['rows']],fresh=fresh,
       costs=json.loads((ROOT/'COST_TOTALS.json').read_text()),
       hypotheses={'finite_step_formula':'DERIVED_STANDARD_RESOLVENT_IDENTITY','local_identifiability_implies_accessibility':'FALSE_BY_ANALYTIC_RATIONAL_VIE_COUNTEREXAMPLE',
       'current_shape_amplitude_ray_opens_failed7_basin':'NOT_SUPPORTED_V1',
       'uniform_multiple_valleys_improve_data_objective':'FINITE_SUPPORT_V2_V3',
       'uniform_multiple_valleys_improve_GT_coverage':'NO_IN_FRESH_3',
       'new_global_convergence_or_probability_solution_set':'NOT_CLAIMED'},
       not_run=['Q6/Q20 proposals (value gate not met)','block Krylov candidate implementation','cross-target uniform Schur cache','new multiple-frequency measurements','fresh reciprocal baseline','global basin/root certification'])
    for name in ('EXPERIMENT_RESULTS.json','RESULTS_SUMMARY.json'):(ROOT/name).write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'compact_endpoints':count,'runs':len(runs),'scored_predictions':len(summary['rows']),'fresh_joint_success':fresh['joint_success_by_method']}))

if __name__=='__main__':main()
