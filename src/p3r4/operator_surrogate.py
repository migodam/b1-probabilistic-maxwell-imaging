"""Truth-free local operator/action supervision from ALREADY-PAID states.

No physical backend is imported. Real <=66 parameter GN solves are separate
from complex Maxwell/current-system RHS. Candidate residuals/fields require an
8-RHS Full state in cold deployment; Full-J labels additionally cost32 RHS.
Predicted errors are empirical scores, never deployment certificates.
"""
import numpy as np
from scipy.special import ndtr
from .core import pack,probit_derivative

BASELINES=('frozenJ','decoder_exact','field_updated','one_secant_Broyden')
ACTION_METHODS=BASELINES+('FullJ',)
ERROR_SCALE=.02
SCALAR_FEATURES=('dimension_over66','log1p_latent_step','material_RMS','linear_material_RMS',
 'decoder_remainder_RMS','new_transition_fraction','anchor_transition_fraction','phase_changed_fraction',
 'prior_energy_delta','Phi_anchor','Phi_new','white_secant_defect_norm','white_linear_action_norm',
 'frozen_J_F','decoder_J_F','field_updated_J_F','Broyden_J_F','D1_material_F',
 'relative_decoder_D_change','relative_total_field_change')
FEATURE_NAMES=tuple([f'residual_new_{i}' for i in range(512)]+[f'residual_anchor_{i}' for i in range(512)]+
 [f'x_new_{i}' for i in range(66)]+[f'delta_x_{i}' for i in range(66)]+list(SCALAR_FEATURES))


def material_matrix(D):
 D=np.asarray(D)
 return np.r_[D.real,D.imag]/np.sqrt(len(D))


def coordinate_derivative(prior,q,x,dimension):
 D=prior.derivative(q)
 return D[:,:64] if dimension==64 else D*probit_derivative(x)[None,:]


def compose_from_transfer(T,E,D,sigma):
 """Exact cached T/total-field/analytic-decoder composition, no inverse solve."""
 T,E,D=map(np.asarray,(T,E,D));d=D.shape[1]
 C=np.einsum('qn,ns->qsn',T,E).reshape(T.shape[0]*E.shape[1],E.shape[0])
 Z=(C@D).reshape(T.shape[0],E.shape[1],d)
 return np.concatenate((Z.real.transpose(1,0,2),Z.imag.transpose(1,0,2)),axis=1).reshape(2*T.shape[0]*E.shape[1],d)*np.sqrt(2)/sigma


def one_secant(J0,h,delta_f):
 h=np.asarray(h);den=float(h@h)
 if den<=1e-24:return J0.copy(),np.zeros(J0.shape[0]),np.zeros_like(h)
 u=np.asarray(delta_f)-J0@h;v=h/den
 return J0+u[:,None]*v[None,:],u,v


def gauss_newton_action(J,r,x,D):
 """H=I+J.TJ, delta=-H^-1(J.Tr+x), material action=D_M delta."""
 J=np.asarray(J,float);r=np.asarray(r,float);x=np.asarray(x,float)
 H=np.eye(J.shape[1])+J.T@J
 delta=np.linalg.solve(H,-(J.T@r+x))  # REAL parameter solve, not a physical complex RHS
 return delta,material_matrix(D)@delta


def ordered_evolution(J0,Jdecoder,Jfield,J1):
 """Exact ordered telescope when T1 was NOT saved.

    Decoder: T0 E0 (D1-D0); field: T0 (E1-E0) D1;
    propagation: (T1-T0) E1 D1 inferred from paid FullJ1-Jfield.
    The last term is teacher-inferred, not recovered T1 or an online feature.
    The alternate T1-based ordering is MISSING_T1 and is never fabricated.
    """
 decoder=np.asarray(Jdecoder)-J0;field=np.asarray(Jfield)-Jdecoder;propagation=np.asarray(J1)-Jfield
 terms=(propagation,field,decoder)
 norms=np.array([float(np.linalg.norm(t)) for t in terms])
 crosses=np.array([float(np.sum(terms[i]*terms[j])) for i,j in ((0,1),(0,2),(1,2))])
 total=J1-J0;sumterms=propagation+field+decoder
 return dict(term_norms=norms,cross_inner_products=crosses,
  telescoping_error=float(np.linalg.norm(sumterms-total)),
  total_defect_F=float(np.linalg.norm(total)),
  squared_norm_identity_error=float(abs(np.linalg.norm(total)**2-(norms@norms+2*crosses.sum()))))


def make_features(*,prior,q0,q1,x0,x1,dimension,D0,D1,E0,E1,J0,Jdecoder,Jfield,Jb,
                  F0,F1,observation,sigma):
 """Only cached anchor + 8-RHS candidate Full state + analytic decoder data.

    No J1, teacher action/error or propagation residual enters the features.
    """
 h=x1-x0;r0=pack(F0-observation)*np.sqrt(2)/sigma;r1=pack(F1-observation)*np.sqrt(2)/sigma
 df=pack(F1-F0)*np.sqrt(2)/sigma
 chi0,a=prior.decode(q0);chi1,b=prior.decode(q1);dchi=chi1-chi0
 phase=lambda f:np.where(f<=-prior.temperature,-1,np.where(f>=prior.temperature,1,0))
 roi=prior.mask;nroi=max(1,int(roi.sum()));M0=material_matrix(D0);M1=material_matrix(D1)
 scalar=np.array([dimension/66,np.log1p(np.linalg.norm(h)),np.linalg.norm(dchi)/24,
  np.linalg.norm(D0@h)/24,np.linalg.norm(dchi-D0@h)/24,
  np.count_nonzero((abs(b['field'])<prior.temperature)&roi)/nroi,
  np.count_nonzero((abs(a['field'])<prior.temperature)&roi)/nroi,
  np.count_nonzero((phase(a['field'])!=phase(b['field']))&roi)/nroi,
  .5*(x1@x1-x0@x0),.5*(r0@r0),.5*(r1@r1),np.linalg.norm(df-J0@h),np.linalg.norm(J0@h),
  np.linalg.norm(J0),np.linalg.norm(Jdecoder),np.linalg.norm(Jfield),np.linalg.norm(Jb),
  np.linalg.norm(M1),np.linalg.norm(D1-D0)/max(np.linalg.norm(D0),1e-30),
  np.linalg.norm(E1-E0)/max(np.linalg.norm(E0),1e-30)])
 xpad=np.zeros(66);hpad=np.zeros(66);xpad[:dimension]=x1;hpad[:dimension]=h
 features=np.r_[r1,r0,xpad,hpad,scalar]
 if features.shape!=(len(FEATURE_NAMES),) or not np.isfinite(features).all():raise ValueError('feature shape/nonfinite')
 return features


def compare_actions(J0,Jdecoder,Jfield,Jb,J1,r,x,D1):
 matrices=(J0,Jdecoder,Jfield,Jb,J1)
 delta=[];action=[]
 for J in matrices:
  a,b=gauss_newton_action(J,r,x,D1);delta.append(a);action.append(b)
 delta=np.stack(delta);action=np.stack(action)
 error=np.linalg.norm(action[:-1]-action[-1],axis=1)
 Jerror=np.array([np.linalg.norm(J-J1)/max(np.linalg.norm(J1),1e-30) for J in matrices[:-1]])
 return dict(delta=delta,material_action=action,material_action_error=error,
  relative_material_action_error=error/max(np.linalg.norm(action[-1]),1e-12),
  relative_J_F_error=Jerror,log_error_target=np.log1p(error/ERROR_SCALE))


def ranker_parameter_count(input_dim=len(FEATURE_NAMES),widths=(256,128,64),output_dim=4):
 sizes=(input_dim,)+tuple(widths)+(output_dim,)
 return sum(a*b+b for a,b in zip(sizes[:-1],sizes[1:]))+sum(2*w for w in widths)


def build_action_error_ranker(input_dim=len(FEATURE_NAMES),widths=(256,128,64)):
 """Lazy CPU PyTorch definition. Defining/preparing a pack performs no training."""
 import torch
 from torch import nn
 layers=[];current=input_dim
 for width in widths:
  layers.extend((nn.Linear(current,width),nn.LayerNorm(width),nn.SiLU()));current=width
 layers.append(nn.Linear(current,len(BASELINES)))
 model=nn.Sequential(*layers)
 if sum(p.numel() for p in model.parameters())>800000:raise ValueError('frozen0.8M parameter cap')
 return model


def empirical_action_errors(log_prediction):
 """Nonnegative empirical estimates, never certified error bounds."""
 return ERROR_SCALE*np.expm1(np.maximum(np.asarray(log_prediction),0))
