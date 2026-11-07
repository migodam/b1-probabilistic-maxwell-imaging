"""Compact final evidence and GN references; no new solves or network calls."""
from pathlib import Path
import json,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1];R4=ROOT.parent/'research_gn_nn_preservation_r4'

def main():
    manifest=json.loads((ROOT/'results/VAL_MANIFEST.json').read_text())
    with np.load(R4/'data/operators.npz')as z:S=np.array(z['S'])
    baseline=[];phases={}
    for phase in ('strict','noiseaware'):
        resultroot=ROOT/'results/noiseaware' if phase=='noiseaware'else ROOT/'results'
        anchors=ROOT/'inputs/noiseaware_anchors'if phase=='noiseaware'else R4/'anchors'
        with np.load(resultroot/'VAL_REFERENCE.npz')as z:
            gt=np.array(z['m_gt']);m=np.array(z['m_gn']);mask=np.array(z['gt_mask'])
        import sys
        sys.path.insert(0,str(ROOT.parents[1]))
        from research.research_gn_search_space_r2.src.score_search import otsu,iou
        selected_ids={r['obs_id']for r in json.loads((ROOT/('audit/noiseaware'if phase=='noiseaware'else'audit')/'CORRECTOR_GROUP_SELECTION.json').read_text())['records']}
        for i,r in enumerate(manifest):
            with np.load(anchors/(r['obs_id']+'.npz'))as z:res=float(np.linalg.norm(S@z['j0']-z['Y'])/np.linalg.norm(z['Y']))
            err=float(np.sqrt(np.mean((m[i]-gt[i])**2)))/r['amplitude'];overlap=iou(m[i]>otsu(m[i]),mask[i]);qualified=res<=(.003 if r['noise'] else .001)
            baseline.append(dict(phase=phase,obs_id=r['obs_id'],parent_id=r['parent_id'],noise=r['noise'],NRMSE=err,IoU=overlap,data_residual=res,qualified=qualified,joint=qualified and err<=.15 and overlap>=.7,corrected_subset=r['obs_id']in selected_ids))
        base_summary=[]
        for subset in ('all156','corrected72'):
            for noise in (0.,.001):
                a=[x for x in baseline if x['phase']==phase and x['noise']==noise and (subset=='all156'or x['corrected_subset'])]
                base_summary.append(dict(scope=subset,noise=noise,cases=len(a),NRMSE=float(np.mean([x['NRMSE']for x in a])),IoU=float(np.mean([x['IoU']for x in a])),qualified=sum(x['qualified']for x in a),joint=sum(x['joint']for x in a)))
        corrected=list(csv.DictReader((resultroot/'CORRECTED_SUMMARY.csv').open()))
        raw=json.loads((resultroot/'DEVELOPMENT_SUMMARY.json').read_text())
        phases[phase]={'baseline_GN':base_summary,'raw_CNN':raw['rows'],'corrected_CNN':corrected}
    with (ROOT/'results/BASELINE_GN_BY_CASE.csv').open('w')as f:w=csv.DictWriter(f,fieldnames=list(baseline[0]));w.writeheader();w.writerows(baseline)
    models=[]
    for folder in ('models','models_noiseaware'):
        for p in sorted((ROOT/folder).glob('*/SUMMARY.json')):
            summary=json.loads(p.read_text());cfg=json.loads((p.parent/'config.json').read_text())
            models.append({'phase':folder,'method':p.parent.name,'epochs':summary['epochs'],'best_epoch':summary['best_epoch'],'parameters':cfg['params'],'physics_counts':summary['physics_counts'],'wall_seconds':summary['wall_seconds']})
    noiseball=json.loads((ROOT/'results/noiseaware/NOISEBALL_RAW_SUMMARY.json').read_text())
    paired=json.loads((ROOT/'results/noiseaware/NOISEBALL_PAIRED_VS_HAAR.json').read_text())
    final={'scope':'DEVELOPMENT TRAIN128/VAL78; no fresh independent TEST','training_observations':256,'validation_observations':156,'independent_validation_mothers':78,'paired_corrector_mothers':36,'phases':phases,'models':models,'noiseball_raw':noiseball,'noiseball_paired_corrected':paired,'scientific_verdict':'Working coupled-loss and fixed-energy CNN prototype; field consistency gain supported in DEV; broad wavelet/fixedpoint/GT superiority HOLD','fresh_TEST':'NOT_RUN','high_contrast':'NOT_RUN','original_FPT_report':'UNAVAILABLE; sourced prior explicit GN fixedpoint equation from R12'}
    (ROOT/'results/RESULTS_SUMMARY.json').write_text(json.dumps(final,indent=2)+'\n')
    print(json.dumps({'model_count':len(models),'baseline':{p:v['baseline_GN']for p,v in phases.items()}},indent=2))

if __name__=='__main__':main()
