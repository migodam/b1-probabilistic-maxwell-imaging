"""Separate additional DEV arm; never overwrites five-arm primary tables."""
from pathlib import Path
import sys,json,csv
import numpy as np
STUDY=Path(__file__).resolve().parents[1];sys.path.insert(0,str(STUDY.parents[1]))
from research.research_gn_search_space_r2.src.score_search import otsu,iou

def main():
    root=STUDY/'results/noiseaware'
    manifest=json.loads((root/'VAL_MANIFEST.json').read_text())
    with np.load(root/'VAL_REFERENCE.npz')as z:
        gt=np.array(z['m_gt']);mask=np.array(z['gt_mask'])
    arm='haar_closed_noiseball_coupled';rows=[]
    with np.load(root/'VAL_PREDICTIONS'/f'{arm}.npz')as z:
        for i,r in enumerate(manifest):
            m=z['m_final'][i];rmse=float(np.sqrt(np.mean((m-gt[i])**2)));nrmse=rmse/r['amplitude']
            overlap=iou(m>otsu(m),mask[i]);qual=bool(z['qualified'][i])
            rows.append(dict(arm=arm,obs_id=r['obs_id'],parent_id=r['parent_id'],noise=r['noise'],NRMSE=nrmse,RMSE=rmse,IoU=overlap,data_residual=float(z['data_residual'][i]),qualified=qual,joint=qual and nrmse<=.15 and overlap>=.7))
    with (root/'NOISEBALL_RAW_BY_CASE.csv').open('w')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=[]
    for noise in (0.,.001):
        a=[r for r in rows if r['noise']==noise]
        summary.append(dict(noise=noise,mothers=len(a),NRMSE=float(np.mean([r['NRMSE']for r in a])),IoU=float(np.mean([r['IoU']for r in a])),data_residual=float(np.mean([r['data_residual']for r in a])),qualified=sum(r['qualified']for r in a),joint=sum(r['joint']for r in a)))
    (root/'NOISEBALL_RAW_SUMMARY.json').write_text(json.dumps({'scope':'additional DEVELOPMENT mechanism arm, notfreshTEST','rows':summary},indent=2)+'\n');print(json.dumps(summary))

if __name__=='__main__':main()
