"""Bounded original-P0 local calibration, never an N3 reconstruction claim.

The only anchors are public TRAIN paired materials. Eight uniform original
Gaussian-shape coordinates vary; all remaining parameters are conditioned on.
No hidden development truth, arbitrary complex pixels, gain, or new prior.
"""
from contextlib import contextmanager
from collections import Counter
from pathlib import Path
import copy,json,time,sys,fcntl,os
import numpy as np
from scipy.special import ndtr,ndtri
from scipy.stats import chi2
from prior_constrained import prior
from prior_constrained.common import pack,ell,descriptor,separation
from .proxy import ROOT,OUT,dump,plain

TAU=float(np.sqrt(chi2.ppf(.95,512)))

class Book:
    def __init__(self):
        self.path=OUT/'COST_LEDGER.jsonl';self.counts=Counter();self.start=time.perf_counter();self.context={}
    def append(self,r):
        with self.path.open('a') as f:f.write(json.dumps(plain(dict(stage='E2',utc=time.time(),**r)))+'\n');f.flush();os.fsync(f.fileno())
    def totals(self):
        a=[json.loads(s) for s in self.path.read_text().splitlines()]
        return dict(states=sum(r.get('amount',0) for r in a if r.get('event')=='STATE'),RHS=sum(r.get('amount',0) for r in a if r.get('event')=='RHS'))
    def reserve(self,event,n,**extra):
        if (OUT/'RUN_CLOSED.json').exists():raise RuntimeError('LXZ stage closed; no automatic new physical authorization')
        with (OUT/'budget.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX);t=self.totals();k='states' if event=='STATE' else 'RHS'
            entries=[json.loads(s) for s in self.path.read_text().splitlines()]
            past_wall=sum(r.get('controlled_wall_seconds',0) for r in entries if r.get('event')=='JOB_END')
            failed=[r for r in entries if r.get('event')=='E2_JOB_FAILURE']
            if failed:
                first=min(r['utc'] for r in entries if r.get('stage')=='E2')
                past_wall+=max(r['utc'] for r in failed)-first
            if time.perf_counter()-self.start+past_wall>3600 or t[k]+n>{'states':96,'RHS':3072}[k]:raise RuntimeError('LXZ user hard physics cap')
            self.append(dict(event=event,amount=n,cumulative=t[k]+n,context=self.context,**extra))
    def reserve_state(self,material,*,resolution,purpose):self.reserve('STATE',1,resolution=resolution,purpose=purpose)
    def reserve_rhs(self,count,role):self.reserve('RHS',int(count),role=role)
    @contextmanager
    def span(self,label,**counts):
        st=time.perf_counter();self.counts.update(counts);status='PASS';error=None
        try:yield
        except BaseException as e:status='FAIL';error=repr(e);raise
        finally:self.append(dict(event='ACTION',label=label,wall_seconds=time.perf_counter()-st,counts=counts,status=status,error=error,context=self.context))

class MaterialPrior:
    """Conditional original smooth-family uniform chart, no support expansion."""
    lo=np.r_[np.full(4,-.13),np.full(2,.05),np.full(2,.4)]
    hi=np.r_[np.full(4,.13),np.full(2,.095),np.full(2,.9)]
    def __init__(self,row,bg):self.row=copy.deepcopy(row);self.bg=bg
    def coordinates(self):
        d=self.row['geometry'];t=np.r_[np.array(d['gaussian_centers'])[:2].ravel(),d['gaussian_widths'][:2],d['gaussian_weights'][:2]]
        return ndtri((t-self.lo)/(self.hi-self.lo))
    def decode(self,u):
        a=self.lo+(self.hi-self.lo)*ndtr(u);r=copy.deepcopy(self.row);d=r['geometry']
        if np.any((a<=self.lo)|(a>=self.hi)):raise FloatingPointError('rounded chart coordinate')
        d['gaussian_centers'][:2]=a[:4].reshape(2,2).tolist();d['gaussian_widths'][:2]=a[4:6].tolist();d['gaussian_weights'][:2]=a[6:].tolist()
        return prior.decode(r,background_id=self.bg),r
    def derivative(self,u):
        h=1e-4;D=[]
        for i in range(8):
            e=np.eye(8)[i]*h;D.append((self.decode(u+e)[0]-self.decode(u-e)[0])/(2*h))
        return np.column_stack(D)
    def analytic_derivative(self,u):
        # Exact derivative within a fixed maximum cell; includes original
        # per-image max normalization and the real-coordinate prior map.
        _,r=self.decode(u);d=r['geometry'];p=prior.geometry_points()-np.array(d['center'])
        c,s=np.cos(d['rotation']),np.sin(d['rotation'])
        xx=(c*p[:,0]+s*p[:,1])/d['scale'];yy=(-s*p[:,0]+c*p[:,1])/(d['scale']*d['axis_ratio'])
        values=[];Dv=np.zeros((len(p),8))
        for i in range(2):
            gx,gy=d['gaussian_centers'][i];w=d['gaussian_weights'][i];a=d['gaussian_widths'][i]
            dx,dy=xx-gx,yy-gy;dist=dx*dx+dy*dy
            g=w*np.exp(-dist/(2*a*a));values.append(g)
            Dv[:,2*i]=g*dx/(a*a);Dv[:,2*i+1]=g*dy/(a*a)
            Dv[:,4+i]=g*dist/a**3;Dv[:,6+i]=g/w
        values=np.sum(values,axis=0);idx=np.argmax(values);peak=values[idx]
        Dv*=((self.hi-self.lo)*np.exp(-.5*u*u)/np.sqrt(2*np.pi))[None,:]
        return 1.2*prior.V*(Dv/peak-values[:,None]*Dv[idx][None,:]/peak**2)

def run():
    dest=OUT/'results/E2_completion';dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE.json').exists():raise RuntimeError('E2 already complete; no silent repeat')
    protocol=json.loads((OUT/'configs/E2_COMPLETION_ERRATUM.json').read_text())
    sys.path.insert(0,str(prior.A30/'src'));from a30.physics import FrozenMaxwellAdapter
    book=Book();allrows=[];allgeom=[];start=time.perf_counter()
    assets=OUT/'delegated/e0_2d'
    for anchor_id in protocol['anchors']:
        meta=json.loads((assets/f'PUBLIC_TRAIN_CALIBRATION_{anchor_id}.json').read_text());bg=anchor_id.split('_')[-1]
        a=np.load(assets/f'PUBLIC_TRAIN_CALIBRATION_{anchor_id}.npz');y=a['observation'];sigma=float(a['sigma'])
        chart=MaterialPrior(meta['generator'],bg);u0=chart.coordinates();x0,row0=chart.decode(u0)
        if np.max(abs(x0-a['chi']))>1e-12:raise RuntimeError('anchor exact decoder mismatch')
        model=FrozenMaxwellAdapter.from_background_npz(prior.INPUTS/'background'/f'{bg}.npz',book,purpose_prefix='LXZ_E2')
        cold=time.perf_counter();book.context=dict(anchor=anchor_id,method='SETUP')
        base=model.full_state(x0,purpose='calibration_anchor_factorization')
        if np.linalg.norm(base.data-a['prediction'])/np.linalg.norm(a['prediction'])>1e-10:raise RuntimeError('cached/actual prediction mismatch')
        transfer=model.receiver_transfer(base);cold=time.perf_counter()-cold
        jac=np.einsum('qp,ps->qsp',transfer,base.field).reshape(256,576)
        D=chart.derivative(u0);C=jac@D
        A=np.column_stack([pack(C[:,i].reshape(32,8))*np.sqrt(2)/sigma for i in range(8)])
        DM=np.r_[D.real,D.imag]/24
        U,s,Vh=np.linalg.svd(DM,full_matrices=False);rank=int(np.sum(s>s[0]*1e-7))
        material_basis=Vh[:rank].T/s[:rank]
        _,sv,Vweak=np.linalg.svd(A@material_basis,full_matrices=False)
        v=material_basis@Vweak[-1];v/=np.sqrt(np.mean(abs(D@v)**2))
        vals,V=np.linalg.eigh(A.T@A);rates=np.maximum(.001,.5/np.sqrt(1+np.maximum(vals,0)))
        B=(V*rates)@V.T;R=(V*np.sqrt(1-rates*rates))@V.T
        e0=ell(base.data,y,sigma);desc0=descriptor(x0,prior.background(bg));rows=[];arrays={};counter=0
        def full(u,method,stage,index,parent=None):
            nonlocal counter
            x,rr=chart.decode(u);book.context=dict(anchor=anchor_id,method=method,stage=stage,index=index)
            st=time.perf_counter();ss=model.full_state(x,purpose=f'{method}_{stage}_{index}');elapsed=time.perf_counter()-st
            support=prior.verify_exact_endpoint(x,rr,background_id=bg)
            en=ell(ss.data,y,sigma);d=float(np.linalg.norm(pack(ss.data-base.data))*np.sqrt(2)/sigma)
            sep=separation(desc0,descriptor(x,prior.background(bg)))
            okay=en<=TAU*TAU/2 and support['prior_status']=='P0_EXACT' and ss.residual<1e-10
            record=dict(anchor=anchor_id,method=method,stage=stage,index=index,id=counter,
                ell=en,delta_ell=en-e0,residual=np.sqrt(2*en),prediction_distance=d,
                material_MAE_to_anchor=float(np.mean(abs(x-x0))),material_RMS_to_anchor=float(np.sqrt(np.mean(abs(x-x0)**2))),
                prior_status=support['prior_status'],compatible=bool(okay),backward_residual=ss.residual,
                wall_seconds=elapsed,posterior_mass=None,candidate_only=method=='continuation',parent=parent,**sep)
            arrays[f'chi_{counter}']=x;arrays[f'prediction_{counter}']=ss.data;arrays[f'u_{counter}']=u
            # Persist every paid endpoint before any later proposal can fail.
            np.savez_compressed(dest/f'{anchor_id}_trial_{counter}.npz',chi=x,prediction=ss.data,u=u,current=ss.current,field=ss.field)
            rows.append(record);dump(dest/f'{anchor_id}_partial_records.json',rows)
            counter+=1;return ss,record
        # Two paid endpoints check actual derivative versus secant. Their values
        # cannot tune the frozen step sizes or select a different target.
        audit=[];h=min(1e-3,.02/np.linalg.norm(v))
        for sign in [-1,1]:audit.append(full(u0+sign*h*v,'derivative_audit','fd',sign)[0])
        fd=(audit[1].data-audit[0].data)/(2*h)
        jv=(C@v).reshape(32,8);derror=float(np.linalg.norm(fd-jv)/max(np.linalg.norm(jv),1e-300))
        for method in ['local','geometry','prior']:
            rng=np.random.default_rng(26101621);u=u0.copy();state=base;en=e0;accepted=[]
            for i in range(protocol['completion_proposals_per_arm']):
                z=rng.normal(size=8)
                if method=='prior':up=z
                elif method=='geometry':up=R@u+B@z
                else:up=np.sqrt(1-.01**2)*u+.01*z
                ss,rec=full(up,method,'proposal',i);ep=rec['ell']
                ok=np.log(rng.random())<min(0,en-ep)
                if ok:u=up;state=ss;en=ep
                rec['MH_accepted']=bool(ok);rec['chain_u']=u.tolist();rec['chain_ell']=en;accepted.append(ok)
            # Repeated rejection states are kept, but this 10-step calibration is
            # explicitly too short for a posterior mass/mixing claim.
        Q=np.eye(8)-np.outer(v,v)/(v@v)
        for i,displacement in enumerate([.003,.006,.012,.024,.048]):
            rawstep=displacement*v;factor=min(1,.5/max(np.linalg.norm(rawstep),1e-300))
            up=u0+factor*rawstep;raw,rec=full(up,'continuation','predictor',i)
            rec['requested_linear_material_displacement']=displacement
            rec['prior_coordinate_trust_factor']=factor
            rec['actual_prior_coordinate_step_norm']=float(np.linalg.norm(up-u0))
            # One normal repair in original parameter coordinates. Keep the
            # selected exploration coordinate; never optimize arbitrary pixels.
            rr=pack(raw.data-y)*np.sqrt(2)/sigma
            AQ=A@Q;delta=-Q@np.linalg.solve(AQ.T@AQ+np.eye(8),AQ.T@rr)
            delta*=min(1,.25/max(np.linalg.norm(delta),1e-300))
            corr,crec=full(up+delta,'continuation','normal_repair',i,parent=rec['id'])
            crec['exploration_coordinate_defect']=float(abs(v@delta)/np.linalg.norm(v))
            crec['raw_displacement_retained']=float(np.mean(abs(corr.chi-x0))/max(np.mean(abs(raw.chi-x0)),1e-300))
        arrays.update(anchor=x0,anchor_prediction=base.data,observation=y,sigma=np.array(sigma),material_Dg=D,whitened_JDg=A,
            material_singular_values=s,response_per_material_singular_values=sv,weak_prior_direction=v)
        np.savez_compressed(dest/f'{anchor_id}.npz',**arrays)
        info=dict(anchor=anchor_id,scope='PUBLIC_TRAIN_KNOWN_MATERIAL_LOCAL_CALIBRATION_ONLY',target='conditional originalP0 smooth chart, eight original uniform parameters; others conditioned fixed',
            cold_setup_seconds=cold,material_chart_rank=rank,latent_dimension=8,material_gauges=8-rank,material_singular_values=s,
            response_per_material_singular_values=sv,prior_response_eigenvalues=vals,derivative_full_secant_relative_error=derror,
            weak_direction_parameter_norm=float(np.linalg.norm(v)),derivative_step=h,
            base_ell=e0,base_residual=np.sqrt(2*e0),noise_sigma=sigma,source_parameters=meta['generator'],
            topology_not_generator_label=True,noise_changed=False,posterior_mass_estimated=False)
        allrows.extend(rows);allgeom.append(info)
        dump(dest/'CANDIDATES.json',dict(rows=allrows,geometry=allgeom))
        print(json.dumps(plain(dict(anchor=anchor_id,chart_rank=rank,derivative_error=derror,compatible=sum(r['compatible'] for r in rows),cost=book.totals()))),flush=True)
    book.append(dict(event='JOB_END',controlled_wall_seconds=time.perf_counter()-start,counts=dict(book.counts)))
    dump(dest/'COMPLETE.json',dict(**book.totals(),wall_seconds=time.perf_counter()-start,GPU=0,full_Mie=0,counts=dict(book.counts)))

if __name__=='__main__':run()
