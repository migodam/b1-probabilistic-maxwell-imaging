"""Disclose seedless synthetic truth and compact predictions after parent scoring closure."""
import json,pathlib,shutil
import numpy as np

P=pathlib.Path(__file__).resolve().parent
S=json.loads((P/'RESULTS_SUMMARY.json').read_text())
F=json.loads((P/S['seal']).read_text())
mapping=json.loads((P/'private_truth/mapping_PRIVATE.json').read_text())
out=P/'truth_after_freeze';out.mkdir(exist_ok=True)
for row in mapping['records']:
    name=row['observation_id']+'.npz'
    shutil.copyfile(P/'private_truth'/name,out/name)
(out/'mapping_PRIVATE.json').write_text(json.dumps({'records':mapping['records'],
    'status':'PUBLIC_SYNTHETIC_TRUTH_DISCLOSED_AFTER_PREDICTION_CLOSURE',
    'original_entropy':'OMITTED','source_seal':S['seal']},indent=2)+'\n')
(out/'README.md').write_text('# Public truth after prediction closure\n\nThese synthetic materials and labels were unavailable to the inverse workers. They are disclosed for offline scoring and reproducibility after all planned dispositions were closed. The historical scorer filename mapping_PRIVATE.json is retained for replay; this copy contains no generator entropy. New blind experiments require new targets and a new freeze.\n')
receipts={r['task_id']:r for r in F['receipts']}
materials=[];valid=[]
for row in S['rows']:
    entry=receipts[row['task_id']]
    result=json.loads((P/entry['path']).read_text()).get('result',{}) if entry.get('path') else {}
    m=result.get('m')
    materials.append(np.asarray(m,float) if row['finite_endpoint'] else np.full(576,np.nan))
    valid.append(row['finite_endpoint'])
(P/'predictions').mkdir(exist_ok=True)
np.savez_compressed(P/'predictions/endpoints.npz',m=np.asarray(materials),finite=np.asarray(valid),
    task_id=np.asarray([r['task_id'] for r in S['rows']]),
    observation_id=np.asarray([r['observation_id'] for r in S['rows']]),
    method=np.asarray([r['method'] for r in S['rows']]),backend=np.asarray([r['backend'] for r in S['rows']]))
