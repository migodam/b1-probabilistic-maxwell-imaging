#!/usr/bin/env python3
"""Accounted native current inversion, with optional exact absorption/reopen."""
from pathlib import Path
import argparse
import json
import os
import sys
sys.dont_write_bytecode=True
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[name]='1'
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from gsda.accounting import BudgetBook,dump
from gsda.fem_model import FEMModel
from gsda.currents import invert,CSIOptions,FixedReference
from gsda.guard import OnlineGuard
from gsda.fwi import FWIOptions


def parse_options(freeze):
    if 'fwi_options' in freeze or 'csi_options' in freeze:
        csi=freeze.get('csi_options',{});prior=freeze.get('fwi_options',{})
    elif 'selected_tv_weight' in freeze:
        csi={};prior=dict(tv_weight=freeze['selected_tv_weight'])
        for name in ('target_chi_scale','target_loss_scale','target_l2_weight','tv_smoothing','wall_spatial_weight'):
            if name in freeze:prior[name]=freeze[name]
        wall=freeze.get('wall_gaussian',{})
        for source,target in (('epsilon_center','wall_eps_center'),('epsilon_sd','wall_eps_sigma'),
                              ('loss_center','wall_loss_center'),('loss_sd','wall_loss_sigma')):
            if source in wall:prior[target]=wall[source]
    else:
        csi=freeze;prior={}
    # Reject malformed option files before constructing a mesh or reserving physics.
    CSIOptions(**csi);FWIOptions(**prior)
    return csi,prior


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--measurement','--observation',dest='measurement',type=Path,required=True)
    ap.add_argument('--output','--out',dest='output',type=Path,required=True)
    ap.add_argument('--options',type=Path)
    ap.add_argument('--job',default='currents')
    ap.add_argument('--absorption',action='store_true')
    ap.add_argument('--risk-controlled',action='store_true',help='B4 frozen empirical screen and actual conflict-triggered REOPEN')
    ap.add_argument('--spectral',action='store_true',help='B3 full-current modal Schur SPD preconditioner with absorption')
    ap.add_argument('--ordinary-schur',action='store_true',help='B5 same full-current modal Schur preconditioner; all material stays active')
    ap.add_argument('--audit-operators',action='store_true',help='Necessary actual-FEM operator/gradient audit only; not inversion')
    args=ap.parse_args()
    if args.ordinary_schur and (args.absorption or args.risk_controlled or args.spectral):
        ap.error('--ordinary-schur is a distinct no-assimilation control and cannot be combined with assimilation flags')
    args.output.mkdir(parents=True,exist_ok=True)
    guard=OnlineGuard(ROOT).install()
    book=BudgetBook(ROOT,job=args.job);status='FAILED'
    try:
        book.check(new_job=True)
        config=json.loads(args.config.read_text());freeze=json.loads(args.options.read_text()) if args.options else {}
        csi,prior=parse_options(freeze)
        if args.spectral or args.risk_controlled or args.ordinary_schur:
            csi={**csi,'spectral_preconditioner':True};CSIOptions(**csi)
        with np.load(args.measurement,allow_pickle=False) as obs:
            if set(obs.files)!={'y'}:raise ValueError('online observation must contain only y')
            y=obs['y'].copy()
        if args.audit_operators:
            # One physical frequency suffices to audit the common quadrature/operator implementation.
            config=dict(config,frequencies=config['frequencies'][:1]);y=y[:1]
        model=FEMModel(config,book,phase='R4_CSI_ORDINARY_SCHUR' if args.ordinary_schur else ('R4_CSI_EMPIRICAL_RISK' if args.risk_controlled else ('R4_CSI_ABSORPTION' if args.absorption or args.spectral else 'R4_CSI')))
        if args.audit_operators:
            from gsda.currents import state_terms
            from gsda.assimilation import absorb_wall,reopen,relative
            ref=FixedReference(model);z=model.initial(0)
            data,js,fields,diag=ref.physical_state(z);rows=[]
            physical=model.predict(z)
            rng=np.random.default_rng(4103002)
            for g,j,field,d in zip(ref.greens,js,fields,data):
                h=(rng.normal(size=j.shape)+1j*rng.normal(size=j.shape))*.01
                gd=g.apply_GD(h)
                physical_current=relative(d-g.ev.direct,g.apply_GS(j))
                state_relative=relative(j,g.chi(z)[:,None]*(g.einc+g.apply_GD(j)))
                jj=j+h;ee=field+gd;alpha=CSIOptions(**csi).state_weight;step=1e-5
                # Complex current derivative includes both GS and GD adjoints.
                chi=g.chi(z);_,_,rr=state_terms(g,z,jj,ee,alpha)
                residual=g.ev.direct+g.apply_GS(jj)-y[g.index]
                gj=g.apply_GSH(residual)/model.sigma**2+alpha*(g.weights[:,None]*rr-g.apply_GDH(np.conj(chi)[:,None]*g.weights[:,None]*rr))
                def objective(jx,ex):
                    rd=g.ev.direct+g.apply_GS(jx)-y[g.index]
                    return .5*float(np.sum(np.abs(rd)**2))/model.sigma**2+state_terms(g,z,jx,ex,alpha)[0]
                fdj=(objective(jj+step*h,ee+step*gd)-objective(jj-step*h,ee-step*gd))/(2*step)
                adjj=float(np.real(np.vdot(gj,h)))
                current_gradient=abs(fdj-adjj)/max(abs(fdj),abs(adjj),1e-12)
                rows.append(dict(frequency=g.p.freq,
                    original_GS0_physical_current_relative=physical_current,state_relative=state_relative,
                    current_gradient_relative=current_gradient,
                    quadrature_points=g.nq,P_shape=g.P.shape,P_nnz=g.P.nnz))
            greens,jr,ea,frozen,absorb=absorb_wall(ref,z,js,fields)
            reopened,ef,reopen_event=reopen(ref,greens,jr,ea,z)
            model_relative=relative(data,physical.data)
            errors=[model_relative]+[v for row in rows for k,v in row.items() if k in ('original_GS0_physical_current_relative','state_relative','current_gradient_relative')]
            result=dict(status='PASS' if max(errors)<1e-4 else 'FAIL',
                scope='one actual same-discrete FEM GS/source consistency, one new complex-current directional adjoint, exact absorption/reopen',
                original_FEM_prediction_relative=model_relative,frequencies=rows,absorption=absorb,reopen=reopen_event,diagnostics=diag)
            dump(args.output/'OPERATOR_AUDIT.json',result)
            if result['status']!='PASS':raise RuntimeError('native current operator audit failed')
        else:
            result=invert(model,y,args.output,csi,prior,args.absorption or args.spectral,args.risk_controlled)
        status=result['status']
    except BaseException as error:
        book.failure('R4_CSI','run_currents',error);dump(args.output/'FAILURE.json',dict(type=type(error).__name__,message=str(error)))
        raise
    finally:
        receipt=book.receipt(status=status,output=str(args.output));dump(args.output/'RECEIPT.json',receipt)
        guard.save(args.output/'FILE_ACCESS_SUMMARY.json')
    print(json.dumps(dict(status=status,output=str(args.output),cost=book.snapshot(pid=os.getpid()))),flush=True)

if __name__=='__main__':main()
