"""Select risk-only fusion on validation; new locked data remain unopened."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,Book,write,save
from p3r2.risk_models import load_predictor
from p3r2.residual_risk import read_anchor,load_residual_predictor
from p3r2.fusion import passive_box,empirical_blend_selection
def main():
    b=Book('D','validation_risk_fusion',max_wall=600,analysis_only=True);status='COMPLETE';error=None
    try:
        a=read_anchor(OUT/'A/anchor_data/validation.npz',True)
        cfg=json.loads((OUT/'configs/TARGET.json').read_text());points=np.array(json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']);roi=np.flatnonzero(np.linalg.norm(points,axis=1)<=cfg['roi_radius'])
        family=[]
        for name,path in [('direct',OUT/'A/training/SELECTED.json'),('scale2',OUT/'A/training32768/SELECTED.json')]:
            if path.exists():
                d=json.loads(path.read_text());p=load_predictor(d['checkpoint']);c=p.predict_full(a['observation']);family.append((name,d['checkpoint'],c))
        residual=json.loads((OUT/'A/residual_training/SELECTED.json').read_text());rp=load_residual_predictor(residual['checkpoint']);family.append(('residual',residual['checkpoint'],rp.predict_full(a)))
        rows=[];pred={}
        for name,path,c in family:
            for projection in ('none','passive_rectangle'):
                v=passive_box(c,cfg,roi) if projection!='none' else c
                selected,curve=empirical_blend_selection(a['chi0'],v,a['chi_GT'])
                key=name+'_'+projection;alpha=selected['alpha'];pred[key]=a['chi0']+alpha*(v-a['chi0'])
                rows.append(dict(name=key,checkpoint=path,projection=projection,alpha=alpha,validation_MSE=selected['complex_MSE'],curve=curve))
        best=min(rows,key=lambda r:(r['validation_MSE'],r['name']))
        write(OUT/'D/fusion/FROZEN.json',dict(selected=best,all_candidates=rows,calibration_parents=a['parent'],GT_scope='explicit supervised validation calibration and model selection, not independent performance evidence',density_or_posterior_claim=False,legal_single_decoder_material=False,freeze_before_new_locked_observations=True,context='trained on48+48contexts; futurestrong96+64covariateshift disclosed'))
        save(OUT/'D/fusion/VALIDATION_PREDICTIONS.npz',parent=a['parent'],**pred)
        print(json.dumps(dict(selected=best['name'],alpha=best['alpha'],validation_MSE=best['validation_MSE'])),flush=True)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)
if __name__=='__main__':main()
