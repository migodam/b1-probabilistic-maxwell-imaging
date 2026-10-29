"""Independent postfreeze scorer. Never imported by the inverse programs."""
import pathlib,json,time,numpy as np
P=pathlib.Path(__file__).resolve().parents[1]
freeze=P/'results/quotient_holdout_FREEZE_RECEIPT.json'
if not freeze.exists():raise RuntimeError('Cannot score unfrozen predictions')
rows=[]
for file in sorted((P/'results').glob('quotient_holdout_QH4_*.json')):
 d=json.loads(file.read_text());case=d.get('case',d.get('task_id'))
 if case is None:raise RuntimeError(list(d))
 truth=np.load(P/'private/quotient_holdout'/f'{case}.npz')['m']
 rr=[]
 for root in d.get('independently_validated',[]):
  if root['status']!='QUALIFIED':continue
  m=np.asarray(root['m']);rr.append({'material_RMS':float(np.linalg.norm(m-truth)/2),'response_relative':root['validation']['response_relative']})
 rows.append({'case':case,'method':d['method'],'qualified_roots':len(rr),'roots':rr,'GT_in_returned_set_RMS_1e5':any(r['material_RMS']<=1e-5 for r in rr)})
if len(rows)!=16:raise RuntimeError(('Expected16 frozen predictions',len(rows)))
out={'status':'SCORED_AFTER_ALL16_PREDICTIONS_FROZEN','scored_unix':time.time(),'freeze_receipt':str(freeze.relative_to(P)),'rows':rows,'scope':'GT used only to assess returned-set containment; no GT-based root selection or parameter changes','GT_values_in_report':False}
(P/'results/quotient_holdout_scores.json').write_text(json.dumps(out,indent=2))
print({'candidate_sets':len(rows),'GT_containment_pass':sum(r['GT_in_returned_set_RMS_1e5'] for r in rows),'max_closest_RMS':max(min(z['material_RMS'] for z in r['roots']) for r in rows)})
