"""Disclosed-target mechanism diagnostic; no hidden material or inverse search."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,importlib.util,json,sys,time
import numpy as np
from scipy.linalg import lu_factor,lu_solve
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent.parent))
from research.research_fbg_r2_adaptive_feedback.src.feedback_fit import fit_feedback,FeedbackOperator
from research.research_fbg_r2_adaptive_feedback.src.schur_current import current_schur,krylov_spaces
from research.research_fbg_r2_adaptive_feedback.src.objective_geometry import directional_geometry

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(ROOT/'results/schur_geometry.json'));a=ap.parse_args()
    old=ROOT.parent/'research_z_alg_r4_diversity'
    with np.load(old/'public_data/operators.npz') as z:G,E,S=z['G'],z['E0'],z['S']
    cases=json.loads((ROOT/'data_contract/stage1_cases.json').read_text())['cases']
    n,p=E.shape;m=np.full(n,.1);L=np.eye(n)-m[:,None]*G
    f=lu_factor(L,check_finite=False);j=lu_solve(f,m[:,None]*E,check_finite=False);Q=E+G@j
    yy,xx=np.meshgrid(np.linspace(0,1,24),np.linspace(0,1,24),indexing='ij')
    hs=np.array([np.ones(n),np.cos(np.pi*xx).ravel(),np.cos(np.pi*yy).ravel(),(np.cos(np.pi*xx)*np.cos(np.pi*yy)).ravel()])
    rhs=(hs[:,:,None]*Q).transpose(1,0,2).reshape(n,-1)
    exact=lu_solve(f,rhs,check_finite=False)
    U,C=krylov_spaces(G,Q,candidate_rank=24)
    t0=time.perf_counter();ret,enriched,detail=current_schur(L,U,C,rhs)
    reduced_time=time.perf_counter()-t0
    T=exact.reshape(n,len(hs),p).transpose(1,0,2)
    Tr=enriched.reshape(n,len(hs),p).transpose(1,0,2)
    Bfit=fit_feedback(G,Q,T,kind='complex',ridge=.05)
    Bschur=fit_feedback(G,Q,Tr,kind='complex',ridge=.05)
    methods={'zero':FeedbackOperator.zero(n),'beta_.2':FeedbackOperator.diagonal(np.full(n,.2+0j)),
             'self':FeedbackOperator.diagonal(-np.diag(G)),'diagonal':Bfit,'Schur_diagonal':Bschur,
             'block':fit_feedback(G,Q,T,kind='block',side=24,ridge=.05)}
    rows=[];beg=time.perf_counter()
    for case in cases:
        with np.load(ROOT/case['public_path']) as z:Y=z['Y']
        R=S@j-Y
        adj=lu_solve(f,S.conj().T@R,trans=2,check_finite=False)
        h=-np.sum(adj.conj()*Q,axis=1).real;h/=np.max(np.abs(h))
        T1=lu_solve(f,h[:,None]*Q,check_finite=False)
        T2=2*lu_solve(f,h[:,None]*(G@T1),check_finite=False)
        for name,fb in methods.items():
            fb,safety=fb.box_certify();B=fb.dense_matrix(G);V=Q+B@j
            geom=directional_geometry(V,R,S,G,B,T1,T2)
            rows.append({'case_id':case['case_id'],'feedback':name,'safety':safety,**geom,
                         'B_diagonal_real':np.diag(B).real.tolist(),'B_diagonal_imag':np.diag(B).imag.tolist()})
    output={'status':'ACTUALLY_RUN_DEVELOPMENT_PUBLIC_DATA','anchor':'uniform .1, same public direction per case acrossfeedbacks',
            'current_core_rank':U.shape[1],'candidate_rank':C.shape[1],'Schur_dimension':C.shape[1],
            'enriched_current_relative_tail':float(np.linalg.norm(exact-enriched)/np.linalg.norm(exact)),
            'retained_current_relative_tail':float(np.linalg.norm(exact-ret)/np.linalg.norm(exact)),
            'fitted_B_relative_difference':float(np.linalg.norm(Bfit.data-Bschur.data)/np.linalg.norm(Bfit.data)),
            'Schur_setup_and_application_wall_s':reduced_time,'diagnostic_wall_s':time.perf_counter()-beg,
            'costs':{'full_LU':1,'source_RHS':p,'exact_current_tangent_RHS':len(hs)*p,
                     'directional_tangent_RHS':2*len(cases)*p,'data_adjoint_RHS':len(cases)*p,
                     'Schur_setup':1,'Schur_RHS':len(hs)*p,'feedback_fits':3},
            'claims':'Galerkin tail measured; exact tangentusedforreference; no A17 runtimeadvantageclaim','rows':rows}
    Path(a.output).write_text(json.dumps(output,indent=2));print(json.dumps({k:output[k] for k in ('enriched_current_relative_tail','fitted_B_relative_difference','Schur_setup_and_application_wall_s')}))
if __name__=='__main__':main()
