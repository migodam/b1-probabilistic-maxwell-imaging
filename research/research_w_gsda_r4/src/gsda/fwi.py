"""Reusable bounded free-material full-wave inversion with spatial TV.

Only public model configuration, measured y, and frozen method configuration
enter this module. World labels, private masks and truth maps are never read.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import asdict, dataclass
import json
import time
import numpy as np
from scipy.optimize import minimize


@dataclass
class FWIOptions:
    # Provisional until METHOD_FREEZE.json is written from independent dev.
    calls_per_stage: tuple=(22,35,150)
    iterations_per_stage: tuple=(15,25,120)
    maxls: int=8
    ftol: float=1e-9
    gtol: float=1e-4
    tv_weight: float=.5
    tv_smoothing: float=.03
    target_chi_scale: float=2.
    target_loss_scale: float=.10
    target_l2_weight: float=.01
    wall_eps_center: float=3.
    wall_eps_sigma: float=1.5
    wall_loss_center: float=.06
    wall_loss_sigma: float=.15
    wall_spatial_weight: float=.05
    frequency_continuation: bool=True
    starts: int=2
    spatial_shapes: object=None

    def __post_init__(self):
        self.calls_per_stage=tuple(int(x) for x in self.calls_per_stage)
        self.iterations_per_stage=tuple(int(x) for x in self.iterations_per_stage)
        if self.starts!=2 or len(self.calls_per_stage)!=3 or min(self.calls_per_stage)<1 or len(self.iterations_per_stage)!=3 or min(self.iterations_per_stage)<1:
            raise ValueError('frozen FWI requires exactly two public starts and three positive stage caps')
        if min(self.tv_smoothing,self.target_chi_scale,self.target_loss_scale,
               self.wall_eps_sigma,self.wall_loss_sigma)<=0:raise ValueError('prior scales must be positive')
        if min(self.tv_weight,self.target_l2_weight,self.wall_spatial_weight)<0:raise ValueError('prior weights must be nonnegative')
        if self.spatial_shapes is not None:
            self.spatial_shapes=tuple(tuple(int(v) for v in shape) for shape in self.spatial_shapes)
            if len(self.spatial_shapes)!=3 or any(len(shape)!=2 or min(shape)<1 for shape in self.spatial_shapes):
                raise ValueError('three positive spatial stage shapes required')


def plain(value):
    if isinstance(value,Path):return str(value)
    if hasattr(value,'tolist'):return value.tolist()
    if hasattr(value,'item'):return value.item()
    raise TypeError(type(value).__name__)


def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,default=plain,indent=2,allow_nan=False)+'\n');temp.replace(path)


def smooth_tv(a,smoothing):
    """Isotropic forward-difference TV, zero outward flux on ROI boundaries."""
    a=np.asarray(a,float);dx=np.zeros_like(a);dy=np.zeros_like(a)
    dx[:,:-1]=a[:,1:]-a[:,:-1];dy[:-1,:]=a[1:,:]-a[:-1,:]
    mag=np.sqrt(dx*dx+dy*dy+smoothing*smoothing)
    px=dx/mag;py=dy/mag;g=np.zeros_like(a)
    g[:,:-1]-=px[:,:-1];g[:,1:]+=px[:,:-1]
    g[:-1,:]-=py[:-1,:];g[1:,:]+=py[:-1,:]
    return float(np.sum(mag-smoothing)),g


def regularization(z,nwall,target_shape,wall_shape,options=None):
    """Shared prior API; returns scalar, full z-gradient and named components."""
    o=FWIOptions() if options is None else options
    if isinstance(o,dict):o=FWIOptions(**o)
    z=np.asarray(z,float);nt=int(np.prod(target_shape));nw=nwall;g=np.zeros_like(z);parts={}
    for name,sl,center,scale in [('wall_eps',slice(0,nw),o.wall_eps_center,o.wall_eps_sigma),
                                ('wall_loss',slice(nw,2*nw),o.wall_loss_center,o.wall_loss_sigma)]:
        q=(z[sl]-center)/scale
        parts[name+'_gaussian']=.5*float(q@q);g[sl]+=q/scale
        a=q.reshape(wall_shape);dx=np.diff(a,axis=1);dy=np.diff(a,axis=0)
        parts[name+'_spatial']=.5*o.wall_spatial_weight*float(np.sum(dx*dx)+np.sum(dy*dy))
        gg=np.zeros_like(a);gg[:,:-1]-=dx;gg[:,1:]+=dx;gg[:-1,:]-=dy;gg[1:,:]+=dy
        g[sl]+=o.wall_spatial_weight*gg.ravel()/scale
    for name,sl,scale in [('target_contrast',slice(2*nw,2*nw+nt),o.target_chi_scale),
                          ('target_loss',slice(2*nw+nt,2*nw+2*nt),o.target_loss_scale)]:
        a=(z[sl]/scale).reshape(target_shape);value,grad=smooth_tv(a,o.tv_smoothing)
        parts[name+'_tv']=o.tv_weight*value
        parts[name+'_l2']=.5*o.target_l2_weight*float(np.sum(a*a))
        g[sl]+=(o.tv_weight*grad.ravel()+o.target_l2_weight*a.ravel())/scale
    return sum(parts.values()),g,parts


def stage_coordinates(model,z,shape):
    """Exact repeated-cell expansion E and Euclidean transpose E.T.

    All twelve wall coordinates remain independent. Coarse target values fill
    public blocks of the same physical material grid; final full shape is the
    identity. Compression is an unweighted block mean for the next-stage seed.
    """
    ny,nx=model.target_shape
    shape=(ny,nx) if shape is None else tuple(shape)
    cy,cx=shape
    if cy>ny or cx>nx:raise ValueError('spatial continuation cannot exceed final material grid')
    iy=np.floor((np.arange(ny)+.5)*cy/ny).astype(int)
    ix=np.floor((np.arange(nx)+.5)*cx/nx).astype(int)
    cells=(iy[:,None]*cx+ix[None,:]).ravel();nb=cy*cx;nw=model.nwall
    mapping=np.r_[np.arange(2*nw),2*nw+cells,2*nw+nb+cells].astype(int)
    ncoord=2*nw+2*nb
    counts=np.bincount(mapping,minlength=ncoord)
    if np.any(counts==0):raise ValueError('empty coarse material coordinate')
    _,first=np.unique(mapping,return_index=True)
    low,high,scale=model.bounds_and_scale()
    q0=np.bincount(mapping,weights=z,minlength=ncoord)/counts
    return mapping,q0,low[first],high[first],scale[first],shape


class CallCap(RuntimeError):
    pass


def invert(model,y,out,options=None):
    o=FWIOptions() if options is None else (FWIOptions(**options) if isinstance(options,dict) else options)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'ONLINE_SEALED.json').exists():raise RuntimeError('refusing to overwrite sealed inversion')
    y=np.asarray(y);nf=len(model.evaluators)
    if y.shape!=(nf,4,16):raise ValueError('online observation shape')
    if o.spatial_shapes is not None and tuple(o.spatial_shapes[-1])!=tuple(model.target_shape):
        raise ValueError('final stage must expose every independent target material cell')
    lo,hi,scale=model.bounds_and_scale();t0=time.time();cpu0=time.process_time()
    stages=[list(range(min(k,nf))) for k in (1,2,nf)] if o.frequency_continuation else [list(range(nf))]*3
    dump(out/'ONLINE_INPUT_CONTRACT.json',dict(observation_keys=['y'],truth_access=False,
        digit_id_access=False,shape_dictionary=False,initialization='two frozen uniform material starts',
        free_target_parameters=2*model.ntarget,free_wall_parameters=2*model.nwall,
        options=asdict(o),config=model.config,stages=stages))
    records=[];totalcalls=0;best=None;last=None;fatal=None
    def save_state(path,state):
        np.savez_compressed(path,z=state['z'],data=state['data'],
            objective=state['objective'],data_objective=state['data_objective'],
            frequency_indices=np.asarray(state['indices']),**model.maps(state['z']))
    try:
        for start in range(o.starts):
            z=model.initial(start);start_records=[]
            for stage,(indices,cap) in enumerate(zip(stages,o.calls_per_stage)):
                calls=0;cache_q=None;cache_result=None;stage_best=None
                shape=None if o.spatial_shapes is None else o.spatial_shapes[stage]
                mapping,q0,stage_lo,stage_hi,stage_scale,stage_shape=stage_coordinates(model,z,shape)
                scaled_lo,scaled_hi=stage_lo/stage_scale,stage_hi/stage_scale
                # Same prior scale relative to expected full-band data objective.
                data_multiplier=nf/len(indices)
                def fun(q):
                    nonlocal calls,totalcalls,cache_q,cache_result,best,last,stage_best
                    if cache_q is not None and np.array_equal(q,cache_q):return cache_result
                    if calls>=cap:raise CallCap('frozen stage objective-call cap')
                    model.book.check();zz=np.clip((q*stage_scale)[mapping],lo,hi)
                    data_value,data_g,pred=model.gradient(zz,y,indices)
                    prior,prior_g,parts=regularization(zz,model.nwall,model.target_shape,model.wall_shape,o)
                    value=data_multiplier*data_value+prior;g=data_multiplier*data_g+prior_g
                    calls+=1;totalcalls+=1
                    coord_g=np.bincount(mapping,weights=g,minlength=len(q))*stage_scale
                    projected=coord_g.copy()
                    projected[((q<=scaled_lo+1e-10)&(coord_g>0))|((q>=scaled_hi-1e-10)&(coord_g<0))]=0.
                    projected_inf=float(np.linalg.norm(projected,np.inf))
                    state=dict(z=zz.copy(),data=pred.data.copy(),objective=value,
                        data_objective=data_value,indices=indices.copy(),start=start,stage=stage,call=calls,
                        scaled_projected_gradient_inf=projected_inf)
                    last=state
                    if stage_best is None or value<stage_best['objective']:stage_best=state
                    comparable=len(indices)==nf
                    if comparable and (best is None or value<best['objective']):
                        best=state;save_state(out/'endpoint.npz',best)
                    save_state(out/'endpoint_last.npz',last)
                    save_state(out/f'start_{start}_stage_{stage}_last.npz',last)
                    if stage_best is state:save_state(out/f'start_{start}_stage_{stage}_best.npz',state)
                    row=dict(time_unix=time.time(),start=start,stage=stage,call=calls,total_calls=totalcalls,
                        objective=value,data_objective=data_value,data_multiplier=data_multiplier,
                        prior=prior,prior_components=parts,frequency_indices=indices,
                        full_objective_comparable=comparable,spatial_shape=stage_shape,optimization_coordinates=len(q),
                        scaled_gradient_norm=float(np.linalg.norm(g*scale)),scaled_projected_gradient_inf=projected_inf,
                        diagnostics=pred.diagnostics)
                    with (out/'TRACE.jsonl').open('a') as handle:handle.write(json.dumps(row,default=plain,allow_nan=False)+'\n')
                    cache_q=q.copy();cache_result=(value,coord_g);return cache_result
                success=False;nit=None;optimizer_status=None;stage_interruption=None
                try:
                    result=minimize(fun,q0/stage_scale,jac=True,method='L-BFGS-B',bounds=list(zip(stage_lo/stage_scale,stage_hi/stage_scale)),
                        options=dict(maxiter=o.iterations_per_stage[stage],maxls=o.maxls,ftol=o.ftol,gtol=o.gtol,maxfun=cap))
                    stop=str(result.message);success=bool(result.success);nit=int(result.nit);optimizer_status=int(result.status)
                    z=stage_best['z'].copy() if stage_best is not None else np.clip((result.x*stage_scale)[mapping],lo,hi)
                except CallCap as error:
                    stop=str(error)
                    if stage_best is not None:z=stage_best['z'].copy()
                except BaseException as error:
                    stop=str(error);stage_interruption=error
                    if stage_best is not None:z=stage_best['z'].copy()
                stage_record=dict(stage=stage,frequency_indices=indices,spatial_shape=stage_shape,optimization_coordinates=len(q0),calls=calls,call_cap=cap,
                    optimizer_success=success,iterations=nit,stop=stop,optimizer_status=optimizer_status,
                    status='OPTIMIZER_TERMINATED_SUCCESS' if success else 'NONCONVERGED',
                    last_scaled_projected_gradient_inf=None if calls==0 or last is None else last.get('scaled_projected_gradient_inf'),
                    best_scaled_projected_gradient_inf=None if stage_best is None else stage_best['scaled_projected_gradient_inf'],
                    last_retained=calls>0,best_retained=stage_best is not None)
                start_records.append(stage_record)
                dump(out/f'START_{start}.json',dict(start=start,stages=start_records))
                if stage_interruption is not None:
                    records.append(dict(start=start,stages=start_records))
                    raise stage_interruption
            records.append(dict(start=start,stages=start_records))
            dump(out/'STARTS.json',records)
    except BaseException as error:
        fatal=dict(type=type(error).__name__,message=str(error))
        if hasattr(model.book,'failure'):model.book.failure(model.phase,'FWI interrupted',error)
        if isinstance(error,(KeyboardInterrupt,SystemExit)):raise
    finally:
        if best is None and last is not None:save_state(out/'endpoint.npz',last)
        result=dict(status='COMPLETE_BOUNDED_ATTEMPT' if fatal is None else 'PARTIAL_INTERRUPTED',
            time_unix=time.time(),case_id=model.config.get('case_id'),method='FULL_FWI_TV_WEIGHTED_ADJOINT',truth_access=False,
            template_used=False,ellipse_used=False,free_target_parameters=2*model.ntarget,
            free_wall_parameters=2*model.nwall,total_objective_calls=totalcalls,
            objective=None if best is None else best['objective'],
            full_objective_comparable=best is not None,
            selected_scaled_projected_gradient_inf=None if best is None else best['scaled_projected_gradient_inf'],
            selection='minimum final full-band measured objective plus frozen TV/spatial prior',
            starts=records,error=fatal,wall_seconds=time.time()-t0,
            process_cpu_seconds=time.process_time()-cpu0,
            retained_last=last is not None,retained_full_objective_best=best is not None,
            convergence_claim='bounded attempts; optimizer stop and nonconvergence retained')
        dump(out/'ONLINE_SEALED.json',result)
    if last is None:raise RuntimeError('no completed physical objective; '+str(fatal))
    return result


def complete_outputs(model,out):
    """Standalone zero-PDE artifacts from an already sealed chosen endpoint."""
    import csv
    out=Path(out)
    if not (out/'ONLINE_SEALED.json').exists():raise RuntimeError('output aliases require online seal')
    with np.load(out/'endpoint.npz',allow_pickle=False) as saved:
        z=saved['z'].copy();data=saved['data'].copy()
        indices=saved['frequency_indices'].tolist() if 'frequency_indices' in saved else list(range(len(model.evaluators)))
    maps=model.maps(z)
    np.save(out/'target_material.npy',np.stack([1+maps['target_contrast'],maps['target_loss']],axis=-1))
    np.save(out/'wall_material.npy',np.stack([maps['wall_eps'],maps['wall_loss']],axis=-1))
    np.savez_compressed(out/'predicted_measurement.npz',data=data,y_pred=data,
        frequencies=np.asarray(model.config['frequencies'])[indices],source_indices=np.arange(4),receiver_indices=np.arange(16))
    rows=[json.loads(line) for line in (out/'TRACE.jsonl').read_text().splitlines() if line]
    fields=['time_unix','start','stage','call','total_calls','objective','data_objective',
            'data_multiplier','prior','frequency_indices','spatial_shape','optimization_coordinates','full_objective_comparable','scaled_gradient_norm','scaled_projected_gradient_inf']
    with (out/'optimization_trace.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader()
        for row in rows:writer.writerow({key:row.get(key) for key in fields})
    dump(out/'reconstruction.json',dict(parameterization='free wall epsilon/loss + free target contrast/loss',
        free_wall_parameters=2*model.nwall,free_target_parameters=2*model.ntarget,
        target_material_format='(ny,nx,2): real permittivity=1+contrast, loss coefficient ell/f',
        wall_material_format='(ny,nx,2): real permittivity, loss coefficient ell/f',
        target_x_edges=maps['target_x_edges'],target_y_edges=maps['target_y_edges'],
        truth_access=False,no_template=True,no_ellipse=True,source='sealed endpoint.npz'))
