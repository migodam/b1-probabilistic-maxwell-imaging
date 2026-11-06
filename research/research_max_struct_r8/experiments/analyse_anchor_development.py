"""Score only the two predeclared, already public validation diagnostics."""
from pathlib import Path
import json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.score_primary import material_metrics

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    manifest=json.loads((ROOT/'configs/warm_anchor_development_manifest.json').read_text())
    records=[]
    for task in manifest['tasks']:
        folder=REPO/task['output_dir']; receipt=json.loads((folder/'RECEIPT.json').read_text())
        if receipt['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError(task['task_id']+' is not sealed')
        records.append((task,folder,receipt))
    # Only the existing validation corpus, never sealed_truth/warm_test.
    with np.load(ROOT/'data_training/validation.npz') as z:truth=z['m'][:8].reshape(8,576).copy()
    cases={c['case_id']:c for c in json.loads((ROOT/'data_public/warm_development/cases.json').read_text())['records']}
    rows=[]
    for task,folder,r in records:
        index=int(task['case_id'].split('_')[-1]);amp=cases[task['case_id']]['condition']['amplitude']
        with np.load(folder/'prediction.npz') as z:m=z['m'].ravel().copy()
        met=material_metrics(m,truth[index],truth[index]>.1+.5*amp,amp)
        rows.append(dict(task_id=task['task_id'],case_id=task['case_id'],arm=task['arm'],initialization=task['initialization'],
            **met,Phi=r['physical_objective'],relative_residual=r['physical_relative_residual'],cpu_seconds=r['cpu_seconds'],
            selected_stage=r['auxiliary_selected_stage'],corrector=r['corrector'],physical_model_cost=r['physical_model_cost'],
            auxiliary_iterations=r['auxiliary_iterations']))
    reference=json.loads((ROOT/'RESULTS/warm_development_analysis/SCORES.json').read_text())['rows']
    reference=[r for r in reference if r['case_id'] in cases and int(r['case_id'].split('_')[-1]) in (0,3)]
    out=ROOT/'RESULTS/warm_anchor_development_analysis';out.mkdir(exist_ok=False)
    (out/'SCORES.json').write_text(json.dumps(dict(scope='TWO_PUBLIC_VALIDATION_OBJECTS; no independent confirmation',rows=rows,existing_reference=reference),indent=2))
    (out/'RECEIPT.json').write_text(json.dumps(dict(status='ANALYSED',cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        LU=0,complex_RHS=0,truth_access='EXISTING_VALIDATION_LABELS_ONLY',primary_truth_access='NONE'),indent=2))
    print(json.dumps([dict(case=r['case_id'],arm=r['arm'],init=r['initialization'],nrmse=r['material_nrmse'],Phi=r['Phi'],cpu=r['cpu_seconds'],stage=r['selected_stage']) for r in rows],indent=2))

if __name__=='__main__':main()
