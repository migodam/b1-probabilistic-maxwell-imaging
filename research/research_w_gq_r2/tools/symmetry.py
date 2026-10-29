from pathlib import Path
from dataclasses import replace
import json,gc,numpy as np
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,WallGrid,Target,S1_TX
from r2.worlds import public_wall,whiteJ,white
ROOT=Path(__file__).resolve().parents[1]
def reflected(d):return np.asarray(d)[[1,0,3,2]][:,::-1]
def main():
    out=ROOT/'results/symmetry';out.mkdir(exist_ok=True);book=BudgetBook()
    w=public_wall(0);eps=w.epsinf.copy();loss=w.ell.copy();xx,yy=np.meshgrid((w.x[:-1]+w.x[1:])/2,(w.y[:-1]+w.y[1:])/2)
    support=w.support;eps[support]=3.2+.25*np.cos(4*xx[support])+.1*(yy[support]-.8)
    # A symmetric pair of air cavities; unknown material perturbations also mirror.
    hole=(abs(xx)<.26)&(yy>.75)&(yy<1);eps[hole]=1.;loss[hole]=0.
    w=WallGrid(w.x,w.y,eps,loss);t=Target([0.,2.1],.12,4.,.07)
    ev=Eval(book,'E1',1.,S1_TX,w.x,w.y);rows=json.loads((out/'RESULTS.json').read_text()) if (out/'RESULTS.json').exists() else []
    basis=[]
    for iy in range(w.shape[0]):
        for ix in range(w.shape[1]//2):
            if w.support[iy,ix] and w.support[iy,-1-ix]:
                col=np.zeros(w.shape);col[iy,ix]=col[iy,-1-ix]=1.;basis.append(col.ravel())
    basis=np.asarray(basis).T
    for case in ['symmetric','asymmetric_wall','differential_gain']:
        if any(r['case']==case for r in rows):continue
        if case=='asymmetric_wall':
            e=w.epsinf.copy();e[w.support&(xx<0)]+=.4;wall=WallGrid(w.x,w.y,e,loss)
        else:wall=w
        state=ev.evaluate(wall,t,label='parity '+case)
        A,B=ev.jacobians(state,basis=basis)
        Ad=A.reshape(4,16,-1);Bx=B[:,0].reshape(4,16)
        oddA=.5*(Ad-Ad[[1,0,3,2],::-1]);evenB=.5*(Bx+reflected(Bx))
        realA=whiteJ(A);realB=white(Bx)
        if case=='differential_gain':
            # Legal shared across-frequency left/right relative amplitude/phase.
            gain=np.array([-1.,1.,-1.,1.])[:,None]*state.data
            realA=np.column_stack((realA,white(gain),white(1j*gain)))
        sol=np.linalg.lstsq(realA,realB,rcond=None)[0];rem=realB-realA@sol
        row=dict(case=case,odd_wall_relative=float(np.linalg.norm(oddA)/np.linalg.norm(Ad)),even_Bx_relative=float(np.linalg.norm(evenB)/np.linalg.norm(Bx)),Bx_norm_white=float(np.linalg.norm(realB)),target_retention=float(np.linalg.norm(rem)**2/np.linalg.norm(realB)**2),singular_values=np.linalg.svd(realA,compute_uv=False).tolist(),rank_scope='machine least-squares only; numerical derivative audit below',instrument='fixed calibrated except explicit relative gain columns')
        fd=[]
        if case=='differential_gain':fd=rows[0]['FD']
        for h in ([] if case=='differential_gain' else [1e-4,.01]):
            dp=ev.evaluate(wall,replace(t,pose=[h,2.1]),label=case+' shiftplus').data.copy()
            dm=ev.evaluate(wall,replace(t,pose=[-h,2.1]),label=case+' shiftminus').data.copy()
            deriv=(dp-dm)/(2*h);fd.append(dict(h=h,relative=float(np.linalg.norm(deriv-Bx)/np.linalg.norm(Bx)),finite_odd_signal=float(np.linalg.norm(.5*(dp-reflected(dp)))/.01)))
            np.savez_compressed(out/f'{case}_{h}.npz',plus=dp,minus=dm,baseline=state.data,analytic_Bx=Bx,A=A,B=B,material=wall.epsinf,loss=wall.ell,x=wall.x,y=wall.y)
        row['FD']=fd;rows.append(row)
        (out/'RESULTS.json').write_text(json.dumps(rows,default=plain,indent=2))
        print(json.dumps(row,default=plain),flush=True)
    print(json.dumps(dict(rows=rows,scope='centered horizontal local task only; target depth/material are not protected',no_rank_certificate=True),default=plain))
if __name__=='__main__':main()
