"""Decisive zero-physics check: new direct local CNN and R2 share a function class."""
from pathlib import Path
import sys,numpy as np,torch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
from p3r4.risk import LocalRisk,load_r2_baseline,feature_image,roi_components,full_material
from p3r2.risk_models import observation_features
OLD=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
def main():
    b=Book('I','cnn_function_class_equivalence',analysis_only=True,max_wall=60)
    status='COMPLETE';error=None
    try:
        torch.set_num_threads(1)
        path=OLD/'A/residual_training/n512/local_geometry_cnn_s33/BEST.pt'
        p=torch.load(path,map_location='cpu',weights_only=False)
        with np.load(OLD/'A/anchor_data/validation.npz',allow_pickle=False) as z:
            online={k:z[k][:8] for k in ['observation','residual','q0','chi0','J']}
        old=load_r2_baseline(path,'cpu');reference=old.predict_full(online)
        model=LocalRisk('local_cnn',p['geometry']);state=model.state_dict()
        for key in state:
            if key.startswith('local.') or key.startswith('head.'):
                state[key]=p['state_dict'][key].clone()
                if key.startswith('head.'):state[key]*=p['residual_scale']
        model.load_state_dict(state);model.eval()
        image=feature_image(observation_features(online['observation']),p['geometry'],online,p['normalization'])
        anchor=roi_components(online['chi0'],p['geometry']['roi'])
        with torch.inference_mode():components=model(torch.tensor(image),torch.tensor(anchor)).numpy()
        actual=full_material(components,p['geometry'])
        error_max=float(np.max(abs(actual-reference)))
        write(OUT/'AUDIT/risk/CNN_CLASS_EQUIVALENCE.json',dict(status='PASS' if error_max<2e-6 else 'FAIL',
            max_complex_error=error_max,parents=8,parameters_new=sum(v.numel() for v in model.parameters()),
            residual_scale=p['residual_scale'],mapping='multiply final head weights and biases by the positive old residual_scale; all three convolutions and eight online input channels identical',
            implication='new local CNN is reoptimized R2 function class, not an independently new architecture or Maxwell geometry mechanism',
            performance_scope='any observed improvement is training/parameterization or data-budget effect until further ablation',
            physics_RHS=0,optimization_updates=0,GT_read=False))
        if error_max>=2e-6:raise AssertionError(error_max)
    except BaseException as e:status='FAILED';error=repr(e);raise
    finally:b.close(status,error)
if __name__=='__main__':main()
