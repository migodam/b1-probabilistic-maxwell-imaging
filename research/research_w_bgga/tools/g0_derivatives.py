"""Paid physical FD/adjoint audit on three predeclared representative states."""
from pathlib import Path
from dataclasses import replace
import gc,json
import numpy as np
from bgga.accounting import BudgetBook,_default
from bgga.fixed_physics import FixedEvaluator
from wslam.forward import Evaluator
from wslam.world import World
from wslam.physics import packj
ROOT=Path(__file__).resolve().parents[1]

def shift(w,name,delta):
    if name.startswith('pose_'):
        p=w.pose.copy();p[0 if name=='pose_x' else 1]+=delta
        return replace(w,pose=p)
    return replace(w,**{name:float(getattr(w,name))+delta})

def main():
    book=BudgetBook();book.authorization('G0')
    out=ROOT/'results/g0_derivatives';out.mkdir(parents=True,exist_ok=True)
    fixed=FixedEvaluator(book)
    reference=Evaluator(book,'G0',h=.17,n=384,grid_refinement=3)
    chosen={'strong_center','right_large','low_loss_high_contrast'}
    cases=[c for c in json.loads((ROOT/'G0_CASES.json').read_text())['qualification'] if c['id'] in chosen]
    names=['epsinf','ell','pose_x','pose_y','contrast','target_loss','radius']
    steps=[1e-3,1e-4,1e-3,1e-3,1e-3,1e-4,1e-3]
    rows=json.loads((out/'RESULTS.json').read_text()) if (out/'RESULTS.json').exists() else []
    for c in cases:
        if any(r['case']==c['id'] for r in rows):continue
        w=World(np.array(c['support'],bool),c['epsinf'],c['ell'],c['pose'],
            radius=c['radius'],contrast=c['contrast'],target_loss=c['target_loss'])
        st=fixed.evaluate(w,label=c['id']);A,B=fixed.derivatives(st)
        J=np.column_stack((A,B[:,:5]));center=st.data.copy();adj_audit=st.derivative_audit
        ref=reference.evaluate(w,jacobian=True);Jref=ref.jacobian.copy();reference.clear()
        fd_columns=[];records=[]
        for name,h0 in zip(names,steps):
            for factor in [1.,.3]:
                h=h0*factor
                def legal(offset):
                    try:
                        ww=shift(w,name,offset)
                        if not .06<=ww.radius<=.16:return None
                        return ww
                    except ValueError:return None
                wm,wp=legal(-h),legal(h)
                if wm is None and wp is None:raise ValueError('No legal finite-difference direction')
                minus=fixed.evaluate(wm,label=f'{c["id"]} {name} -{h}').data.copy() if wm is not None else None
                plus=fixed.evaluate(wp,label=f'{c["id"]} {name} +{h}').data.copy() if wp is not None else None
                one_sided=wm is None or wp is None
                fd=(plus-center)/h if wm is None else (center-minus)/h if wp is None else (plus-minus)/(2*h)
                jcol=J[:,names.index(name)]
                records.append(dict(parameter=name,h=h,one_sided=one_sided,stencil=('forward' if wm is None else 'backward' if wp is None else 'central'),
                    fd_relative=float(np.linalg.norm(fd.ravel()-jcol)/np.linalg.norm(jcol)),
                    absolute_whitened_column_error=float(np.linalg.norm(packj((fd.ravel()-jcol)[:,None])))))
                if factor==.3:fd_columns.append(fd.ravel())
        fd=np.column_stack(fd_columns)
        row=dict(case=c['id'],duality=adj_audit,finite_differences=records,
            fixed_vs_fitted_J_relative=float(np.linalg.norm(J-Jref)/np.linalg.norm(Jref)),
            empirical_operator_error=float(np.linalg.norm(packj(J-Jref),2)),
            empirical_FD_operator_error=float(np.linalg.norm(packj(J-fd),2)),
            error_scope='measured discrepancies, not uniform derivative bounds',
            names=names)
        rows.append(row)
        np.savez_compressed(out/(c['id']+'.npz'),J=J,J_reference=Jref,J_FD=fd,A=A,B=B,data=center)
        (out/'RESULTS.json').write_text(json.dumps(rows,default=_default,indent=2)+'\n')
        print(json.dumps(row,default=_default),flush=True)
        del st;gc.collect()

if __name__=='__main__':main()
