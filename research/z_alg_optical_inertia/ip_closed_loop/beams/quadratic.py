"""Old-conditioned quadratic forms and source-index proposals (float proposals)."""
import numpy as np
from scipy.linalg import lu_factor, lu_solve, eigh

from beam_common import herm


class ResidualForm:
    """Use q=J-chi_center E, so E=E_center a+Tq and Y=Y_center a+Dq."""
    def __init__(self, widths, field, response, t, d):
        self.m, self.p = field.shape[1], d.shape[0]
        self.field, self.response, self.t, self.d = field,response,t,d
        self.indices = np.arange(len(widths))
        self.eye=np.eye(len(widths))
        wt=widths[:,None]*t
        we=widths[:,None]*field
        self.base_a = herm(self.eye-wt.conj().T@wt)
        self.base_n = -(wt.conj().T @ we)
        self.base_c = -(we.conj().T @ we)
        self.readout_right=d.conj().T@np.column_stack((response,-np.eye(self.p)))
        self.readout_ext=herm(np.column_stack((response,-np.eye(self.p))).conj().T @
                              np.column_stack((response,-np.eye(self.p))))
        self.row_norm=np.linalg.norm(t,axis=1)

    def prepare(self, cells, mu, *, state):
        n = len(self.eye)
        u=np.column_stack((self.eye[:,cells], self.t[cells,:].conj().T))
        base=self.base_a*mu
        right=self.base_n*mu
        constant=self.base_c*mu
        if state:
            base=base+self.d.conj().T@self.d
            right=np.column_stack((right,np.zeros((n,self.p),complex)))+self.readout_right
            constant=np.pad(constant,((0,self.p),(0,self.p)))+self.readout_ext
        inv=lu_factor(base,check_finite=False)
        dim=right.shape[1]
        solved=lu_solve(inv,np.column_stack((right,u)),check_finite=False)
        return (cells, mu, right, constant, u, solved[:,:dim], solved[:,dim:])

    def update(self, prepared, i, ratio):
        cells,mu,right,const,u,kr,ku=prepared
        width=len(cells)
        cell=cells[i]
        u_cols=[i,i+width]
        ur=u[:,u_cols]
        ud=ku[:,u_cols]
        t=-.5j*mu*ratio
        mult=np.array([[0,t],[np.conj(t),0]],complex)
        little=herm(np.linalg.inv(mult)+ur.conj().T@ud)
        v=eigh(little,eigvals_only=True,check_finite=False)
        nu=int(np.sum(v>0))-1
        e=np.zeros(right.shape[1],complex)
        e[:self.m]=t*self.field[cell]
        z=right+ur[:,0,None]*e[None,:]
        zk=kr+ud[:,0,None]*e[None,:]
        w=zk-ud@np.linalg.solve(little,ur.conj().T@zk)
        return nu,herm(const-z.conj().T@w)


def old_border(q,m,y,old_radius):
    old=np.r_[np.eye(m)[:,0],y]
    compound=q@old
    scalar=float((np.vdot(old,compound)-old_radius**2).real)
    nu=int(scalar>0)
    return scalar,nu,herm(q-np.outer(compound, compound.conj())/scalar)


def envelope(q,m,vector):
    r=q[m:,m:]
    rm=np.linalg.solve(r,q[m:,:m])
    source=herm(q[:m,:m]-q[:m,m:]@rm)
    excess=-float(np.vdot(vector,source@vector).real)
    return r, -rm@vector, excess,source


def quad_search(form, cells, y, radius_old, radius_new, gram):
    """Return frozen candidates, with source/receiver diagnostics for every trial."""
    powers=(.015625,.125,1.,8.,64.,512.,4096.)
    ratios=(1.8, 2.05, 2.2, 2.5, 3., 4., 8.)
    physical=form.prepare(cells,1.,state=False)
    ps={(i,u):form.update(physical,i,u/form.row_norm[cell])
        for i,cell in enumerate(cells) for u in ratios}
    all_trials=[]
    candidates=[]
    for mu in powers:
        system=form.prepare(cells,mu,state=True)
        for i,cell in enumerate(cells):
            for u in ratios:
                ratio=u/form.row_norm[cell]
                nphys,source=ps[i,u]
                nu,q=form.update(system,i,ratio)
                record=dict(cell=int(cell),state_negative=nu,physical_negative=nphys,mu=mu,ratio=ratio)
                if nu!=1 or nphys!=1:
                    all_trials.append(record)
                    continue
                old,k,conditioned=old_border(q,form.m,y,radius_old)
                r=conditioned[form.m:,form.m:]
                evr=eigh(r,eigvals_only=True,check_finite=False)
                evals,evectors=eigh(source[1:,1:],gram[1:,1:],check_finite=False)
                positive=max(0., float(evals[-1])*mu)
                min_needed=radius_new**2/positive if positive>0 else None
                record.update(old_positive=k, source_positive=int(sum(evals>0)),old_pivot=float(old),
                              generalized_source_top=float(evals[-1]),energy_threshold_for_one_source=min_needed,
                              minimum_receive_eig=float(evr[0]))
                all_trials.append(record)
                if k==nu and evr[0]>1.e-10:
                    score=float(evr[0])
                    # The source-only eigenbeam is an explicit first-class proposal.
                    source_mode=np.r_[0.,evectors[:,-1]]
                    candidates.append((score,record,conditioned,source_mode,source))
    candidates.sort(key=lambda item:item[0],reverse=True)
    return candidates,all_trials

