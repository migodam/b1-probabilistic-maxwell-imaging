"""Finite reciprocal probe, all passive loads, original terminated line array."""
from pathlib import Path
from time import perf_counter, process_time
import json, sys
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent
PRIOR=HERE.parents[1]/'z_alg_optical_inertia'/'blind_uq_mimo_theory'
sys.path.insert(0,str(PRIOR/'ports'))
from model import DiskGreen, network_terms, point
NEW_POSITION=np.array([.12,1.12])
ACTIVE_PORTS=np.r_[np.arange(16),24]
R=50.


class Circuit:
    def __init__(self,side,band,active,position=NEW_POSITION):
        self.side,self.band=side,band
        self.disk=DiskGreen(side,band)
        opfile=PRIOR/'ports/public/operators'/(
            f'n{side}_f400_a08.npz' if band==400 else f'n{side}_f{band}_a08_c048.npz')
        with np.load(opfile) as d:
            self.symbol=d['spectrum']; self.pad=tuple(d['pad'])
            old_b,old_c,z0,old_pos,disk_pos,old_gain=(d[k] for k in
                    ('E0_bare','S_bare','Z_AA','positions','centers','gain'))
            pulse=float(d['pulse'])
        self.ids=np.asarray(active)
        yy,xx=np.unravel_index(self.ids,(side,side))
        lag=np.fft.ifft2(self.symbol)
        self.g=lag[(yy[:,None]-yy)%self.pad[0],(xx[:,None]-xx)%self.pad[1]]
        feed=np.sqrt(48.)*point(np.linalg.norm(disk_pos[self.ids]-position,axis=1),self.disk.k)
        self.b=np.column_stack((old_b[self.ids],feed))
        self.c=np.vstack((old_c[:,self.ids],pulse*feed))
        coupling=-1j*old_gain*np.sqrt(48.)*point(np.linalg.norm(old_pos-position,axis=1),self.disk.k)
        self.z0=np.zeros((25,25),dtype=complex)
        self.z0[:24,:24]=z0
        self.z0[:24,24]=self.z0[24,:24]=coupling
        self.z0[24,24]=20.+0j
        self.pulse=pulse
        self.count={'lu':0,'rhs':0,'jacobian':0}

    def solve(self,x,jac=False):
        a=np.eye(len(x),dtype=complex)-x[:,None]*self.g
        lu=lu_factor(a,check_finite=False)
        j=lu_solve(lu,x[:,None]*self.b,check_finite=False)
        z=self.z0-1j*(self.c@j)
        inv=np.linalg.solve(z+R*np.eye(25),np.eye(25))
        total=np.eye(25)-2*R*inv
        self.count['lu']+=1;self.count['rhs']+=25
        if not jac: return total
        e=self.b+self.g@j
        c=lu_solve(lu,self.c.T,trans=1,check_finite=False).T
        self.count['rhs']+=25;self.count['jacobian']+=1
        return total,(inv@c),e@inv

    def raw_with_residual(self,x):
        a=np.eye(len(x),dtype=complex)-x[:,None]*self.g
        lu=lu_factor(a,check_finite=False)
        return lu_solve(lu,x[:,None]*self.b,check_finite=False), lu


def at_load(s,k,ports=tuple(range(16))):
    """k is the power-wave reflection at the added, undriven 25th terminal."""
    ports=np.asarray(ports)
    return s[np.ix_(ports,ports)]+k*np.outer(s[ports,24],s[24,ports])/(1-k*s[24,24])


def summarize_pair(sp,sm,r):
    from itertools import combinations
    tx=np.array([(0,*t) for t in combinations(range(1,16),7)])
    rx=np.array([[i for i in range(16) if i not in t] for t in tx])
    k=np.exp(2j*np.pi*np.arange(256)/256)
    result={}
    for key,gamma in [('open',1.),('matched',0.),('short',-1.)]:
        d=at_load(sp,gamma)-at_load(sm,gamma)
        cut=d[rx[:,:,None],tx[:,None,:]]
        fro=np.linalg.norm(cut,axis=(1,2))/np.sqrt(8)
        result[key]=dict(min_cut=float(fro.min()/(2*r)),max_cut=float(fro.max()/(2*r)),
                         median_cut=float(np.median(fro)/(2*r)),
                         A=float(np.linalg.norm(d[8:16,:8])/np.sqrt(8)/(2*r)),
                         B=float(np.linalg.norm(d[1::2,::2])/np.sqrt(8)/(2*r)),
                         arbitrary_wave_spectral_ratio=float(np.linalg.norm(d,2)/(2*r)))
    d=np.array([at_load(sp,v)-at_load(sm,v) for v in k])
    result['all_loads_grid_sup_transfer_over_2r']=float(np.linalg.norm(d,ord=2,axis=(1,2)).max()/(2*r))
    return result


def material_task(side):
    ax=(np.arange(side)+.5)*2/side-1
    yy,xx=np.meshgrid(ax,ax,indexing='ij')
    rho=np.hypot(xx-.125,yy-.0625).ravel()
    ix=np.flatnonzero(rho<=.50)
    edge=np.linspace(.15,.5,13)
    group=np.searchsorted(edge,rho[ix],side='left')
    return ix,group


def refine(side):
    """Fit the complete matched active-plus-probe scattering matrix, no target truth or data."""
    from scipy.optimize import least_squares
    start=perf_counter(),process_time()
    ix,group=material_task(side)
    wave=Circuit(side,400,ix)
    states=0
    upper=4.5
    with np.load(PRIOR/'ports/public/cases'/f'n{side}_c01.npz') as d:r=float(d['radius'])
    last={}
    def val(y):
        nonlocal states
        if 'y' in last and np.array_equal(y,last['y']): return
        xp,xm=np.zeros(len(ix)),np.zeros(len(ix))
        xp[group==0]=3.5
        for i,v in enumerate(y,1):
            xp[group==i]=upper/2+v
            xm[group==i]=upper/2-v
        p,lp,rp=wave.solve(xp,True)
        m,lm,rm=wave.solve(xm,True)
        dif=(p-m)[np.ix_(ACTIVE_PORTS,ACTIVE_PORTS)]
        j=np.stack([(-2j*R*((lp[:,group==i]@rp[group==i])+(lm[:,group==i]@rm[group==i])))[np.ix_(ACTIVE_PORTS,ACTIVE_PORTS)]
                    for i in range(1,len(y)+1)],axis=-1).reshape(dif.size,-1)
        last.update(y=y.copy(),diff=np.r_[dif.real.ravel(),dif.imag.ravel()],jac=np.vstack((j.real,j.imag)),p=p,m=m,xp=xp,xm=xm)
        states+=1
    def fun(y):val(y);return last['diff']
    def jac(y):val(y);return last['jac']
    # Repeat at native sampling, using the physical shell amplitudes only as a start.
    if side>32:
        with np.load(HERE/'adversary_32.npz') as d:initial=d['amplitudes'].copy()
    else: initial=np.zeros(12)
    fun(np.zeros(12))
    uncomp=summarize_pair(last['p'],last['m'],r)
    opt=least_squares(fun,initial,jac=jac,bounds=(-upper/2,upper/2),method='trf',
                      gtol=1e-11,ftol=1e-10,xtol=1e-9,max_nfev=160)
    val(opt.x)
    result=summarize_pair(last['p'],last['m'],r)
    xplus,xminus=np.zeros(side**2),np.zeros(side**2)
    xplus[ix]=last['xp'];xminus[ix]=last['xm']
    p,m=last['p'],last['m']
    cost400=dict(wave.count)
    other=Circuit(side,460,ix)
    p460=other.solve(last['xp']); m460=other.solve(last['xm'])
    np.savez_compressed(HERE/f'adversary_{side}.npz',ix=ix,group=group,amplitudes=opt.x,
        material_plus=xplus,material_minus=xminus,plus400=p,minus400=m,plus460=p460,minus460=m460,
        stored_probe_excitation=wave.b[:,-1],stored_probe_readout=wave.c[-1,:],augmented_empty_Z=wave.z0)
    ext=xplus-xminus
    ext[ix[group==0]]=0
    born=[]
    for x,s in [(last['xp'],p),(last['xm'],m)]:
        op0=np.linalg.inv(wave.z0[:24,:24]+R*np.eye(24))
        born0=(-2j*R*op0@(wave.c[:24]@(x[:,None]*wave.b[:,:24]))@op0)[:16,:16]
        signal=(at_load(s,1)-(np.eye(24)-2*R*op0)[:16,:16])[8:16,:8]
        born.append(float(np.linalg.norm(born0[8:16,:8]-signal)/np.linalg.norm(signal)))
    out=dict(side=side,reference_noise_radius_sqrtW=r,core_cells=int(np.sum(group==0)),
        core_area_m2=float(np.sum(group==0)*(2/side)**2),physical_exterior_L1=float(sum(abs(ext))*(2/side)**2),
        physical_exterior_L2=float(np.linalg.norm(ext)*2/side),physical_exterior_Linf=float(max(abs(ext))),
        n_unknowns=side**2,external_group_counts=np.bincount(group,minlength=13).tolist(),
        nfev=opt.nfev,models_evaluated=states,stop=opt.message,optimality=opt.optimality,
        compensation_400=result,uncompensated_400=uncomp,untrained_460=summarize_pair(p460,m460,r),
        born_relative_error_A=born,cost_400=cost400,cost_460=other.count,
        wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1])
    (HERE/f'adversary_{side}.json').write_text(json.dumps(out,indent=2)+'\n')
    print('grid',side,'nfev',opt.nfev,'A/B',result['open']['A'],result['open']['B'],'all_load_upper_grid',result['all_loads_grid_sup_transfer_over_2r'],'cost',out['wall_s'],flush=True)


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('sides',nargs='+',type=int);args=ap.parse_args()
    with threadpool_limits(1):
        for side in args.sides:refine(side)
