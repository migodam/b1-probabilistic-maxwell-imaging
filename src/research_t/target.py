"""One frozen actual A30 backend; no forward computation on import."""
from pathlib import Path
import json,sys,os,csv
import numpy as np
from scipy.stats import chi2
from .budget import dump

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'research/research_t_20261015'
A30=Path(os.environ.get('A30_ROOT',str(ROOT.parents[1]/'A30/a30-maxwell-material-manifold')))
RELEASE=A30/'inputs/a29_release'
CASES=[f'p{p}_{b}_a1_snr20' for p in [250065,250069,250073,250077] for b in ['weak','ring']]
TAU=float(np.sqrt(chi2.ppf(.95,512)))

def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).reshape(-1)
def state_arrays(state):return dict(chi=np.array(state.chi),current=np.array(state.current),field=np.array(state.field),prediction=np.array(state.data),backward_residual=np.array(state.residual))
def save_state(path,state):np.savez_compressed(path,**state_arrays(state))

class Targets:
    def __init__(self,book):
        self.book=book;self.models={}
        src=str(A30/'src')
        if src not in sys.path:sys.path.insert(0,src)
        from a30.physics import FrozenMaxwellAdapter
        self.Adapter=FrozenMaxwellAdapter
    def case(self,case):
        if case not in CASES:raise ValueError('outside frozen exposedcase list')
        bg='weak' if '_weak_' in case else 'ring'
        if bg not in self.models:
            allowed=[RELEASE/'inputs/a28/results/a28/exact_states'/f'{c}.npz' for c in CASES if f'_{bg}_' in c]
            self.models[bg]=self.Adapter.from_background_npz(RELEASE/'inputs/a25/background'/f'{bg}.npz',self.book,allowed_cache_paths=allowed,purpose_prefix='RESEARCH_T')
        with np.load(RELEASE/'results/a29/pools'/f'{case}.npz',allow_pickle=False) as z:
            # Online fields are deliberately enumerated. No truth/teacher loads.
            obs=z['observation'].copy();sigma=float(z['sigma']);n3=z['n3'].copy()
        anchor=self.models[bg].cached_anchor(n3,RELEASE/'inputs/a28/results/a28/exact_states'/f'{case}.npz')
        return self.models[bg],anchor,obs,sigma

def validate(state,observation,sigma):
    x=state.chi;res=float(np.linalg.norm(pack(state.data-observation)*np.sqrt(2)/sigma))
    con=float(np.linalg.norm(state.current-x[:,None]*state.field)/max(np.linalg.norm(state.current),1e-300))
    legal=bool(np.all(x.real>=-.5-1e-12) and np.all(x.real<=4+1e-12) and np.all(x.imag>=-1e-12) and np.all(x.imag<=.5+1e-12))
    okay=bool(res<=TAU and legal and con<1e-10 and state.residual<1e-10)
    return dict(status='FULL_VALIDATED_NUMERICAL' if okay else 'REJECTED',residual_norm=res,tau=TAU,material_valid=legal,constitutive_defect=con,backward_residual=float(state.residual),output_error_upper_bound=None,numerical_bound_status='UNKNOWN_NUMERICAL_ERROR',shared_material=True,nuisance=dict(gain=1,known=True,shared_sources=8))
