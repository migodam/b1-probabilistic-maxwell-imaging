"""Paid, explicitly scoped audit of R5's sufficient chart-radius bound.

The norm/residual formula is an exact-arithmetic bound. Its float64 evaluation
is NOT a verified-roundoff or continuum certificate. Centers use observations
only, selected by fixed IDs from the retained v1 optimization; no truth access.
"""
import os
for key in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import svdvals
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from research.p3l_r5_parallel_ab.shared.physics_adapter import Target,SIGMA,save,write

OUT=Path(__file__).resolve().parent/'curvature_bounds'
HIST=OUT.parents[1]/'stage_a_response_coordinates/v1_failed_preserved/physical'

def run():
    OUT.mkdir(exist_ok=True)
    if (OUT/'RESULTS.json').exists():
        raise RuntimeError('Already completed; no implicit paid rerun')
    cfg=dict(version=1,centers=['S0_P0_N0','S0_P1_N0'],mode_index=0,
             R=3.,noise_multipliers=[1.,.5,.25],forward_states=2,
             forward_rhs=16,inverse_rhs=1152,jvp_rhs=32,
             norm_scope='float64 audit of exact-arithmetic global decoder/Neumann bound',
             truth_used=False)
    write(OUT/'FROZEN.json',cfg)
    target=Target('A','parent_curvature_bounds',output_dir=OUT/'physics')
    before=target.book.totals();rows=[]
    for ident in cfg['centers']:
        mode=json.loads((HIST/(ident+'_MODES.json')).read_text())['modes'][0]
        x=np.zeros(66);x[:2]=mode['z']
        p=target.evaluate(x,purpose='A5_observation_only_center_'+ident,retain_reference=True)
        V=np.eye(66)[:,:2]
        J=SIGMA*target.jacobian(p,V,purpose='A5_direction_frame_'+ident)
        _,sing,vt=np.linalg.svd(J,full_matrices=False)
        W=V@vt[0];kappa=float(sing[0])
        if not kappa>0:raise RuntimeError('No arbitrary spectral floor')
        R0=target.inverse_action(p,np.eye(576,dtype=complex),purpose='A5_inverse_norm_'+ident,
             persist_path=OUT/(ident+'_INVERSE.npz'))
        with target.book.compute('derivative','A5_global_norm_algebra_'+ident):
            model=target.backend();prior=target.prior
            g=float(svdvals(model.G)[0]);s=float(svdvals(model.S)[0])
            defect=float(np.linalg.norm(np.eye(576)-p.state.matrix@R0,'fro'))
            if defect>=1:raise RuntimeError('Inverse defect cannot give stability bound')
            rnorm=float(svdvals(R0)[0])
            beta=rnorm/(1-defect)
            beta_fro=float(np.linalg.norm(R0,'fro')/(1-defect))
            # Global bounds along the fixed, unit, scalar strong direction W.
            coeff=(prior.A@W[:64])/prior.temperature
            knorm=float(np.max(abs(coeff[prior.mask])))
            dr=float(np.diff(prior.real_interval)[0]);di=float(np.diff(prior.imag_interval)[0])
            delta_max=float(abs(complex(prior.real_interval[1],prior.imag_interval[1])))
            delta1=float(np.hypot(dr*W[64],di*W[65])/np.sqrt(2*np.pi))
            delta2=float(np.hypot(dr*W[64]**2,di*W[65]**2)/np.sqrt(2*np.pi*np.e))
            h1=15/16;h2=2.5/np.sqrt(3)
            D1=delta1+delta_max*h1*knorm
            D2=delta2+2*delta1*h1*knorm+delta_max*h2*knorm**2
            E0=float(np.linalg.norm(p.state.field,'fro'))
            def L_of(rho):
                perturb=D1*rho;eta=beta*g*perturb
                if eta>=1:return float('inf')
                b=beta/(1-eta)
                ef=(1+g*b*perturb)*E0
                return np.sqrt(2)*s*b*ef*(D2+2*g*b*D1**2)
            safe=.49/(beta*g*D1)
            rho=brentq(lambda r:L_of(r)*r-kappa/2,0,safe,xtol=1e-15)
            L=L_of(rho)
            for mul in cfg['noise_multipliers']:
                sigma=SIGMA*mul;R=cfg['R']
                rows.append(dict(center=ident,noise_multiplier=mul,kappa=kappa,
                    G_norm=g,S_norm=s,R0_spectral_norm=rnorm,inverse_defect_F=defect,
                    beta_formula=beta,beta_Frobenius_formula=beta_fro,
                    decoder_strong_D1_global_bound=D1,decoder_strong_D2_global_bound=D2,
                    sufficient_rho=rho,L_on_rho_ball=L,
                    curvature_radius_ratio=L*rho/kappa,
                    noise_radius_ratio=2*sigma*R/(kappa*rho),
                    L_sigma_R_over_kappa_squared=L*sigma*R/kappa**2,
                    max_sufficient_response_R=kappa*rho/(2*sigma),
                    requested_R=R,local_sufficient_conditions_at_R3=bool(2*sigma*R<=kappa*rho),
                    level='EXACT_ARITHMETIC_FORMULA_EVALUATED_FLOAT64_NOT_CERTIFIED'))
            save(OUT/(ident+'_FRAME.npz'),x=x,V=V,W=W,J=J,chi=p.chi,prediction=p.prediction,
                 mode_source=np.array(str(HIST/(ident+'_MODES.json'))))
        del p,R0
    after=target.book.totals()
    result=dict(status='COMPLETED',config=cfg,rows=rows,cost_before=before,cost_after=after,
                own_intended_rhs=1200,
                cost_note='Shared A totals can include interleaved worker actions; job receipts identify own cost.',
                interpretation='Failure of this sufficient radius is not failure of the exact inverse or nonexistence of posterior uncertainty.',
                limitations=['computed matrix norms have no directed-rounding enclosure',
                  'centers need not have zero selected residual; this audits radius necessities only',
                  'not a full posterior coverage or weak-marginal certificate'])
    write(OUT/'RESULTS.json',result)
    print(json.dumps(dict(status=result['status'],rows=len(rows),conditions=[r['local_sufficient_conditions_at_R3'] for r in rows])))

if __name__=='__main__':run()
