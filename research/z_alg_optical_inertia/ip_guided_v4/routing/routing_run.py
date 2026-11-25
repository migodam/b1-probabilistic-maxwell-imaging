#!/usr/bin/env python3
"""Public-data Born / full-wave routing on the inherited Maxwell operators."""
from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import minimize
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
RESEARCH = STUDY.parent
ROOT = RESEARCH.parent
BORN = STUDY / 'physical_revision/born'
LOSS = STUDY / 'ip_unified/numerics/loss'
SRAP = STUDY / 'validation/srap'
RING = STUDY / 'validation/max_struct'
NATIVE_SOLVER = RESEARCH / 'z_alg_xudong_cie/minpeak_native'
sys.path[:0] = [str(BORN), str(NATIVE_SOLVER)]
import extinction
from flint import ctx
import search as native
from curve_inputs import full_ring_inputs


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    entry = importlib.util.module_from_spec(spec)
    sys.modules[name] = entry
    spec.loader.exec_module(entry)
    return entry


loss = module('routing_loss', LOSS / 'run.py')


def read(path):
    return json.loads(path.read_text())


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False, default=lambda x: x.item()) + '\n')


def clock():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])


@dataclass
class Observation:
    id: str
    object_id: str
    cohort: str
    geom: str
    source: str
    y: np.ndarray
    radius: float
    noise: float
    port_status: str = 'withheld simulated ports scored offline'


@dataclass
class Geometry:
    s: np.ndarray
    b: np.ndarray
    g: np.ndarray
    side: int
    width: float
    frequency: float
    dense: bool = False
    tangent: float = 0.
    lo: float = 0.
    hi: float = 4.


def native_suite(include_2pct=False, only=None):
    spec, s, b, original = full_ring_inputs()
    geometry = Geometry(s, b, spec, 96, 2., 400e6)
    base = STUDY / 'submission/curve/observations'
    for r in original:
        if (only is None or r.id in only or r.object_id in only) and (r.noise_pct == 10 or (include_2pct and r.noise_pct == 2)):
            yield geometry, Observation(r.id, r.object_id, 'optfresh', 'native_10x20',
                                        str((base / (r.id + '.npz')).relative_to(ROOT)),
                                        r.y, r.padded_radius, r.noise_pct)
    with np.load(SRAP / 'public_observations.npz', allow_pickle=False) as data:
        online = {q['id']: q for q in read(SRAP / 'online.json')['objects']}
        for ident, y in zip(data['ids'], data['Y'], strict=True):
            object_id = str(ident.split('__')[0])
            if only is None or ident in only or object_id in only:
                yield geometry, Observation(str(ident), object_id, 'srap22', 'native_10x20',
                                            str((SRAP/'public_observations.npz').relative_to(ROOT)),
                                            y, online[ident]['padded_radius'], 10.)


def other_suite(which):
    if which == 'weak':
        geom = next(native_suite())[0]
        base = BORN/'controls'
        with np.load(base/'arrays.npz', allow_pickle=False) as data:
            for row in read(base/'results.json')['observations']:
                name = row['observation']
                obj = name.split('__')[0].replace('_born', '_vie')
                yield geom, Observation(name,obj,'v3_weak_controls','native_10x20',
                                        str((base/'arrays.npz').relative_to(ROOT)),
                                        data[name+'__Y'],row['radius'],100*row['nominal_noise_fraction'],
                                        'no frozen holdout; offset full-wave test only for VIE controls')
    elif which == 'ring':
        with np.load(RING/'public_inputs.npz', allow_pickle=False) as data:
            width = float(np.ptp(data['centers'][:,0])*24/23)
            geom = Geometry(data['S'],data['E0'],data['G'],24,width,None,True,lo=.02,hi=20.)
            for name, y, radius in zip(data['cases'],data['Y'],data['delta'],strict=True):
                yield geom, Observation(str(name)+'__24to24-12tx-48rx',str(name),'r14_ring7','ring_12x48',
                                        str((RING/'public_inputs.npz').relative_to(ROOT)),y,float(radius),1.,
                                        'no frozen separately sampled ports')
    elif which == 'loss':
        for i,width in enumerate((.6,1.25)):
            base = loss.finite_geometry(width)
            file = LOSS / f'finite_{i}_inputs.npz'
            with np.load(file, allow_pickle=False) as data:
                for row in read(LOSS/f'finite_{i}_verified.json')['inputs']:
                    name = row['id']
                    obj,ts,_ = name.split('__')
                    if obj.endswith('_low'):
                        continue
                    t = float(ts.replace('fixed_t',''))
                    geom = Geometry(base.s,base.b,base.dense,12,width,400e6,True,t)
                    yield geom, Observation(name,obj,'loss4',f'loss_{i}_t{t:g}',str(file.relative_to(ROOT)),
                                            data[name],row['rho'],2.)
    elif which == 'mild_mid':
        f = HERE/'mild_mid_public.npz'
        with np.load(f, allow_pickle=False) as data:
            geom = Geometry(data['S'],data['B'],data['G'],12,.6,400e6,True)
            for name, y, r in zip(data['ids'],data['Y'],data['radii'], strict=True):
                yield geom, Observation(str(name),str(name).split('__')[0],'new_mild_mid','mild_mid_10x20',
                                        f.name,y,float(r),2.)
    else:
        raise ValueError(which)


def offset_geometry(geom, receive_radius):
    k = 2*np.pi*geom.frequency/3e8
    x = (np.arange(geom.side)+.5)*(geom.width/geom.side)-geom.width/2
    yy,xx = np.meshgrid(x,x,indexing='ij')
    points = np.c_[xx.ravel(),yy.ravel()]

    def directions(num):
        angle = 2*np.pi*(np.arange(num)+.5)/num
        return np.c_[np.cos(angle),np.sin(angle)]

    d = np.linalg.norm(receive_radius*directions(len(geom.s))[:,None]-points[None],axis=2)
    kr = k*geom.width/geom.side/np.sqrt(np.pi)
    s = .5j*np.pi*kr*jv(1,kr)*hankel1(0,k*d)
    b = np.exp(1j*k*(points@directions(geom.b.shape[1]).T))
    return Geometry(s,b,geom.g,geom.side,geom.width,geom.frequency,geom.dense,geom.tangent,geom.lo,geom.hi)


class DenseModel:
    def __init__(self, geo, y):
        self.geo, self.y = geo, y
        self.scale = np.linalg.norm(y)**2
        self.count = Counter()

    def full(self, m, gradient=False):
        geo = self.geo
        complex_factor = 1+1j*geo.tangent
        x = complex_factor*m
        n = len(x)
        a = -x[:,None]*geo.g
        a.flat[::n+1] += 1
        factors = lu_factor(a,check_finite=False)
        self.count['native_full_wave_calls'] += 1
        self.count['state_rhs'] += geo.b.shape[1]
        j = lu_solve(factors,x[:,None]*geo.b,check_finite=False)
        p = geo.s@j
        r = p-self.y
        fr = float(np.linalg.norm(r))
        if gradient:
            self.count['adjoint_rhs'] += geo.b.shape[1]
            adj = lu_solve(factors,geo.s.conj().T@r,trans=2,check_finite=False)
            grd = np.sum((complex_factor*(geo.b+geo.g@j)*adj.conj()).real,axis=1)/self.scale
            return .5*fr*fr/self.scale, grd, fr, j, p
        return fr,p,j


def full_audit(geo, y, m):
    started = clock()
    if geo.dense:
        model = DenseModel(geo,y)
        fr,ports,j = model.full(m)
        defect = j - ((1+1j*geo.tangent)*m)[:,None]*(geo.b+geo.g@j)
    else:
        model = native.Model(geo.g,geo.s,geo.b,y)
        rtol, native.RTOL = native.RTOL, 1.e-10
        try:
            fr,ports,j = model.full(m,gradient=False)
        finally:
            native.RTOL = rtol
        defect = j-m[:,None]*(geo.b+model.green(j))
    return dict(residual=float(fr), peak=float(np.max(m)), legal=bool(np.all(np.isfinite(m)) and
                np.min(m)>=geo.lo and np.max(m)<=geo.hi),
                state_defect=float(np.linalg.norm(defect)), counts=dict(model.count), **elapsed(started)),ports,j


def fit_dense(geo, obs, initial, calls):
    started = clock()
    model = DenseModel(geo,obs.y)
    best = {}

    def fun(m):
        f, grad, residual, current, ports = model.full(m,gradient=True)
        if residual < best.get('residual',np.inf):
            best.update(residual=residual,m=m.copy())
        if residual < .998*obs.radius:
            raise native.Feasible
        return f,grad

    try:
        res = minimize(fun,initial,jac=True,method='L-BFGS-B',bounds=[(geo.lo,geo.hi)]*len(initial),
                       options={'maxiter':calls,'maxfun':calls,'gtol':1.e-10,'ftol':2.e-12,'maxls':25,'maxcor':10})
        stop = str(res.message)
    except native.Feasible:
        stop = 'reached numerical interior of noise ball'
    return best['m'],dict(stop=stop,counts=dict(model.count),**elapsed(started))


def fit_iterative(geo, obs, initial, calls):
    if geo.dense:
        return fit_dense(geo,obs,initial,calls)
    model = native.Model(geo.g,geo.s,geo.b,obs.y)
    state, stat = native.box_fit(model,initial,np.arange(24**2),24,geo.hi,
                                 .998*obs.radius,calls,120.,True)
    return state['m'],{**stat,'counts':dict(model.count)}


def born_setup(geo):
    t = clock()
    basis = np.einsum('ri,is->rsi',geo.s,geo.b).reshape(-1,len(geo.b))
    if not geo.dense:
        index = native.groups(24)
        basis = np.array([np.bincount(index,weights=r.real,minlength=24**2) +
                          1j*np.bincount(index,weights=r.imag,minlength=24**2) for r in basis])
    a = basis*(1+1j*geo.tangent)
    h = (a.conj().T@a).real
    lip = float(np.linalg.eigvalsh(h)[-1])
    return (a,h,lip),elapsed(t)


def born_fit(geo, obs, setup):
    started = clock()
    a,h,lip = setup
    y = obs.y.ravel()
    linear = (a.conj().T@y).real
    ridge = .001*lip
    z = x = np.full(len(h),geo.lo)
    alpha = 1.
    for _ in range(120):
        nxt = np.clip(z-(h@z-linear+ridge*z)/(lip+ridge),geo.lo,geo.hi)
        new_alpha = .5*(1+np.sqrt(1+4*alpha**2))
        z = nxt+(alpha-1)/new_alpha*(nxt-x)
        x,alpha = nxt,new_alpha
    ports = (a@x).reshape(obs.y.shape)
    m = x if geo.dense else x[native.groups(24)]
    return m,ports,dict(residual=float(np.linalg.norm(ports-obs.y)), iterations=120, **elapsed(started))


def certificate_setup(geo, ident):
    started = clock()
    b = geo.b*(1+1j*geo.tangent)
    if ident == 'native_10x20':
        with np.load(BORN/'fixed_port_weights.npz',allow_pickle=False) as old:
            z = old['w_10_20']
        meta = {'proposal':'unchanged v3 fixed weights'}
    else:
        meta,z = extinction.weights(geo.s,b)
        meta['proposal'] = 'same deterministic v3 geometry rule on declared operator'
    if z is None:
        return None, None, {**meta,'valid':False,**elapsed(started)}
    zw,n,margin = extinction.geometry_check(geo.s,b,z)
    output = dict(**meta,valid=bool(margin<0),max_cell_upper=float(np.nextafter(float(margin.upper()),np.inf)),**elapsed(started))
    return (zw,n) if output['valid'] else None,z,output


def cert_readout(ob, cert):
    t = clock()
    if cert is None:
        return dict(status='no verified fixed polar direction',lower=None,excludes=False,**elapsed(t))
    v = extinction.readout(*cert,ob.y)
    return dict(status='strict Born cone excluded' if v>ob.radius else 'fixed lower bound inconclusive',
                lower=v,excludes=v>ob.radius,**elapsed(t))


def mild_data():
    base = loss.finite_geometry(.6)
    geom = Geometry(base.s,base.b,base.dense,12,.6,400e6,True)
    x,y = loss.positions(12,.6).T
    profiles = [(.002, (y>-.14)&(y<.17-1.4*abs(x+.05))),
                (.03, ((abs(x-.07)<.035)&(abs(y)<.18))|((abs(y+.06)<.035)&(abs(x)<.18))),
                (.25, ((x+.06)**2+(y-.06)**2<.145**2)&((x-.04)**2+(y-.14)**2>.081**2)),
                (.8, ((abs(x+.10)<.09)&(abs(y+.065)<.115))|((abs(x-.115)<.065)&(abs(y-.115)<.055)))]
    ids, fields, radii, materials = [], [], [], []
    port = np.zeros((20,10),complex)
    for i,(amp,area) in enumerate(profiles):
        name = f'guided_mildmid_{i+1:02d}__pw360_10x20__pair2pct'
        m = amp*area
        _,clean,_ = DenseModel(geom,port).full(m)
        rng = np.random.default_rng(np.random.SeedSequence([20261125088,i]))
        phase = rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape)
        dy = .02*np.linalg.norm(clean)*phase/np.linalg.norm(phase)
        y = clean+dy
        ids.append(name)
        fields.append(y)
        radii.append(np.linalg.norm(y-clean)+1.e-6*np.linalg.norm(y))
        materials.append(m)
    np.savez_compressed(HERE/'mild_mid_public.npz',S=geom.s,B=geom.b,G=geom.g,ids=ids,Y=fields,radii=radii)
    np.savez_compressed(HERE/'mild_mid_offline_reference.npz',ids=ids,materials=materials)


def run_one(geo, obs, linear, cert, out, calls):
    tr = cert_readout(obs,cert)
    born, born_ports, t_born = born_fit(geo,obs,linear)
    paths = [('born_only',born),('full_zero',np.full(len(geo.b),geo.lo)),('born_then_full',born)]
    if sum(map(ord,obs.id))%2 == 0:
        # The physical audit remains immediately after each arm. Rotate the two nonlinear starts.
        paths[1:] = paths[:0:-1]
    steps, arrays = {}, {'Y':obs.y,'born_prediction':born_ports}
    for tag,start in paths:
        whole = clock()
        try:
            m, fit = (start,None) if tag=='born_only' else fit_iterative(geo,obs,start,calls)
            audit,p,j = full_audit(geo,obs.y,m)
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as error:
            steps[tag] = dict(failure=str(error), **elapsed(whole))
            continue
        arrays.update({tag+'_material':m,tag+'_ports':p,tag+'_current':j})
        steps[tag] = dict(fit=fit,audit=audit,**elapsed(whole))
    born_acceptable = (t_born['residual']<=obs.radius and
                       steps['born_only'].get('audit',{}).get('residual',np.inf)<=obs.radius)
    if tr['excludes']:
        selected = 'full_zero'
        tail = steps[selected]['wall_seconds']
        tail_cpu = steps[selected]['cpu_seconds']
    elif born_acceptable:
        selected = 'born_only'
        tail = t_born['wall_seconds']+steps['born_only']['wall_seconds']
        tail_cpu = t_born['cpu_seconds']+steps['born_only']['cpu_seconds']
    else:
        selected = 'born_then_full'
        tail = sum(r['wall_seconds'] for r in (t_born,steps['born_only'],steps[selected]))
        tail_cpu = sum(r['cpu_seconds'] for r in (t_born,steps['born_only'],steps[selected]))
    np.savez_compressed(out/(obs.id+'.npz'),**arrays)
    result = dict(id=obs.id,object_id=obs.object_id,cohort=obs.cohort,setup_key=obs.geom,source=obs.source,
                  forward_grid=geo.side,inversion_grid=geo.side if geo.dense else 24,width_m=geo.width,
                  frequency_Hz=geo.frequency,noise_pct=obs.noise,radius=obs.radius,
                  field_norm=float(np.linalg.norm(obs.y)),n_tx=obs.y.shape[1],n_rx=obs.y.shape[0],
                  lo=geo.lo,hi=geo.hi,tangent=geo.tangent,holdout=obs.port_status,born_certificate=tr,
                  born_linear_fit=t_born,arms=steps,
                  route=dict(selected=selected,valid_born_and_physical_point=born_acceptable,
                             charged_wall_seconds=tail+tr['wall_seconds'],charged_cpu_seconds=tail_cpu+tr['cpu_seconds']))
    write(out/(obs.id+'.json'),result)
    print(obs.id, selected, [round(steps[k].get('audit',{}).get('residual',np.inf)/obs.radius,4)
                             for k in ['born_only','full_zero','born_then_full']],flush=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--suite',choices=['native','weak','ring','loss','mild_mid'],required=True)
    parser.add_argument('--calls',type=int,default=16)
    parser.add_argument('--out',type=Path,default=HERE/'online')
    parser.add_argument('--objects',nargs='+')
    parser.add_argument('--include-2pct',action='store_true')
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    if args.suite == 'mild_mid' and not (HERE/'mild_mid_public.npz').exists():
        mild_data()
    cache, metadata = {}, {}
    ctx.prec,ctx.threads = 80,1
    started = clock()
    cases = native_suite(args.include_2pct,args.objects) if args.suite == 'native' else other_suite(args.suite)
    done = []
    for geo,obs in cases:
        if obs.geom not in cache:
            cert,z,proof = certificate_setup(geo,obs.geom)
            linear,time = born_setup(geo)
            metadata[obs.geom] = dict(born_geometry=proof,born_fit_setup=time)
            cache[obs.geom] = linear,cert
            if z is not None:
                np.savez_compressed(args.out/(args.suite+'__setup_'+obs.geom+'.npz'),z=z)
        run_one(geo,obs,*cache[obs.geom],args.out,args.calls)
        done.append(obs.id)
    write(args.out/(args.suite+'__seal.json'),dict(ids=done,setup=metadata,calls=args.calls,**elapsed(started)))


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
