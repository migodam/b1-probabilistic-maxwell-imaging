"""Offline exposed-VAL score, only after all public predictions are saved."""
from pathlib import Path
import sys,json,csv,argparse
import numpy as np
STUDY=Path(__file__).resolve().parents[1];sys.path.insert(0,str(STUDY.parents[1]))
from research.research_gn_search_space_r2.src.score_search import otsu,iou

ARMS=['pixel_closed_coupled','haar_closed_coupled','db2_closed_coupled','haar_closed_no_fullfield','haar_free_coupled']
R4=STUDY.parent/'research_gn_nn_preservation_r4'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--anchor-mode',choices=['strict','noiseaware'],default='strict');args=parser.parse_args()
    root=STUDY/'results/noiseaware' if args.anchor_mode=='noiseaware' else STUDY/'results'
    anchor_root=STUDY/'inputs/noiseaware_anchors' if args.anchor_mode=='noiseaware' else R4/'anchors'
    root.mkdir(parents=True,exist_ok=True)
    records=[r for r in json.loads((R4/'data/manifests/evaluation.json').read_text())['records']if r['split']=='val']
    ids=[r['obs_id']for r in records]
    for arm in ARMS:
        p=root/'VAL_PREDICTIONS'/f'{arm}.npz'
        if not p.exists():raise RuntimeError('Public predictions incomplete; do not score partial method set')
    gt=[];gn=[];masks=[]
    for r in records:
        with np.load(R4/'data'/r['truth_path'])as z:
            gt.append(np.array(z['m']).reshape(576));masks.append(np.array(z['alpha']).reshape(576)>=.5)
        with np.load(anchor_root/(r['obs_id']+'.npz'))as z:gn.append(np.array(z['m0']).reshape(576))
    gt=np.array(gt);gn=np.array(gn);masks=np.array(masks)
    ref=dict(obs_id=np.array(ids),parent_id=np.array([r['parent_id']for r in records]),noise=np.array([r['noise']for r in records]),amplitude=np.array([r['amplitude']for r in records]),m_gt=gt,m_gn=gn,gt_mask=masks)
    np.savez_compressed(root/'VAL_REFERENCE.npz',**ref)
    (root/'VAL_MANIFEST.json').write_text(json.dumps(records,indent=2)+'\n')
    rows=[]
    for arm in ARMS:
        with np.load(root/'VAL_PREDICTIONS'/f'{arm}.npz')as z:
            if list(z['obs_id'])!=ids:raise ValueError('Prediction order must match declared public VAL records')
            for i,r in enumerate(records):
                m=z['m_final'][i];err=float(np.sqrt(np.mean((m-gt[i])**2)));rel=float(z['data_residual'][i]);qual=bool(z['qualified'][i])
                threshold=otsu(m);overlap=float(iou(m>threshold,masks[i]))
                rows.append(dict(arm=arm,obs_id=r['obs_id'],parent_id=r['parent_id'],split='DEVELOPMENT_VAL',noise=r['noise'],amplitude=r['amplitude'],family=r['family'],digit=r.get('digit'),font=r.get('font'),material_NRMSE=err/r['amplitude'],material_RMSE=err,IoU=overlap,data_residual=rel,qualified=qual,joint=qual and err/r['amplitude']<=.15 and overlap>=.7,inference_seconds_amortized=float(z['inference_seconds_amortized'][i]),forward_wall_seconds=float(z['inference_seconds_amortized'][i]),seen_image=False,seen_mother=False,seen_class=True,seen_font=True))
    out=root/'RESULTS_BY_CASE.csv'
    with out.open('w')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=[]
    for arm in ARMS:
        for noise in (0.,.001):
            a=[r for r in rows if r['arm']==arm and r['noise']==noise]
            summary.append(dict(arm=arm,noise=noise,observations=len(a),mothers=len({r['parent_id']for r in a}),NRMSE=float(np.mean([r['material_NRMSE']for r in a])),RMSE=float(np.mean([r['material_RMSE']for r in a])),IoU=float(np.mean([r['IoU']for r in a])),data_residual=float(np.mean([r['data_residual']for r in a])),qualified=sum(r['qualified']for r in a),joint=sum(r['joint']for r in a)))
    (root/'DEVELOPMENT_SUMMARY.json').write_text(json.dumps({'scope':'DEVELOPMENT VAL used for checkpoints; not fresh TEST','anchor_mode':args.anchor_mode,'rows':summary,'new_TEST':'NOT_RUN'},indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
