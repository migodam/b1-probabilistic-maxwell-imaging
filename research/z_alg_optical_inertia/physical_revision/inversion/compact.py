"""Prospective electrically compact full-wave transition experiment."""

import json
import sys

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

from fit import native, run
from gate import HERE, STUDY, Observation, clocks, elapsed, write

sys.path.insert(0, str(STUDY))
from core.certificate import ball, certify, frobenius_squared, safe_float
from core.inertia import Prepared

NATIVE = 12
WIDTH = .125
K = 2*np.pi*400e6/3e8
LOW, HIGH = .5, 2.
GAMMAS = (.45, .60, .75, .95)
WEIGHT_SCALES = (.03, .10, .30, 1., 3., 10.)
STAGES = ((6, 20), (9, 30), (12, 44))
OUT = HERE/'compact_synthetic'


def unit(angles):
    return np.stack((np.cos(angles), np.sin(angles)), axis=-1)


def acquisition(held=False):
    x = (np.arange(NATIVE)+.5)*WIDTH/NATIVE-WIDTH/2
    yy, xx = np.meshgrid(x,x,indexing='ij')
    points = np.stack((xx.ravel(), yy.ravel()),axis=-1)
    tx = (np.arange(10)+.5*held)*2*np.pi/10
    rx = (np.arange(20)+.5*held)*2*np.pi/20
    radii = WIDTH / NATIVE / np.sqrt(np.pi)
    dist = np.linalg.norm(points[:,None]-points[None,:],axis=2)
    outgoing = lambda d: .5j*np.pi*K*radii*jv(1,K*radii)*hankel1(0,K*d)
    g = outgoing(np.where(dist>0, dist, 1.))
    np.fill_diagonal(g,.5j*np.pi*K*radii*hankel1(1,K*radii)-1)
    s = outgoing(np.linalg.norm(.375*unit(rx)[:,None]-points[None,:],axis=2))
    b = np.exp(1j*K*(points@unit(tx).T))
    return (g,s,b), points


def material(points):
    x,y = (2*points/WIDTH).T
    rr = (np.abs(x-.06)/.95)**10+(np.abs(y+.035)/.92)**10
    envelope = np.clip((1.1-rr)/.45,0.,1.)
    raise_lobe = np.exp(-((x+.28)/.40)**2-((y-.24)/.37)**2)
    lower_lobe = np.exp(-((x-.31)/.28)**2-((y+.22)/.34)**2)
    return envelope*(.78+.1*raise_lobe-.06*lower_lobe)


class DenseModel(native.Model):
    def __init__(self, *args):
        super().__init__(*args)
        self.forward_m = None

    def green(self, z, adjoint=False):
        self.count['dense_G_products'] += z.shape[1]
        return (self.spectrum.conj().T if adjoint else self.spectrum)@z

    def linear_solve(self, material, rhs, adjoint, seed):
        if not np.any(material):
            return rhs.copy()
        if self.forward_m is None or not np.array_equal(material,self.forward_m):
            self.factor = lu_factor(np.eye(len(material))-material[:,None]*self.spectrum,check_finite=False)
            self.forward_m = material.copy()
            self.count['dense_state_factorizations'] += 1
        self.count['adjoint_rhs' if adjoint else 'state_rhs'] += rhs.shape[1]
        return lu_solve(self.factor,rhs,trans=2 if adjoint else 0,check_finite=False)


def enclosure(geometry, m):
    """Full stored dense model, enclosed with an independent inverse-defect bound."""
    start = clocks()
    g,s,b = geometry
    n = len(m)
    a0 = np.eye(n)-m[:,None]*g
    fac = lu_factor(a0,check_finite=False)
    raw_inverse = ball(lu_solve(fac, np.eye(n),check_finite=False))
    ball_a = acb_mat(n,n,[ (int(i==j)-arb(float(m[i]))*v)
                           for i,row in enumerate(g) for j,v in enumerate(ball(row[None,:]).entries())])
    ball_f = acb_mat(*b.shape,[arb(float(m[i]))*z
                              for i,row in enumerate(b) for z in ball(row[None,:]).entries()])
    identity = ball(np.eye(n))
    defect = frobenius_squared(identity-ball_a*raw_inverse).sqrt().upper()
    if defect>=1:
        raise ArithmeticError('The inverse residual did not prove invertibility.')
    inv_bound = frobenius_squared(raw_inverse).sqrt().upper()/(1-defect)
    currents = ball(lu_solve(fac,m[:,None]*b,check_finite=False))
    field_defect = frobenius_squared(ball_f-ball_a*currents).sqrt().upper()
    sensors = ball(s)
    remainder = frobenius_squared(sensors).sqrt().upper()*inv_bound*field_defect
    ports = sensors*currents
    return ports, remainder, dict(inverse_defect_upper=safe_float(defect,upper=True),
                                  port_error_upper=safe_float(remainder,upper=True), **elapsed(start))


def norm(z):
    return frobenius_squared(z).sqrt()


def observe():
    start = clocks()
    public, private = OUT/'public', OUT/'offline_truth'
    public.mkdir(parents=True,exist_ok=True)
    private.mkdir(parents=True,exist_ok=True)
    geometry, points = acquisition()
    held, _ = acquisition(held=True)
    truth = material(points)
    p, rem, audit = enclosure(geometry,truth)
    held_ports, held_rem, held_audit = enclosure(held,truth)
    clean = np.array([complex(z) for z in p.entries()]).reshape(20,10)
    pristine_held = np.array([complex(z) for z in held_ports.entries()]).reshape(20,10)
    clean_legal_p, legal_rem, _ = enclosure(geometry,np.full(len(truth),LOW))
    gap = np.linalg.norm(clean-np.array([complex(z) for z in clean_legal_p.entries()]).reshape(20,10))
    rng = np.random.default_rng(2026112403)
    base = rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape)
    base /= np.linalg.norm(base)
    records = []
    for gamma in GAMMAS:
        sample = gamma*gap
        y = clean+sample*base
        field = ball(y)
        radius = safe_float(norm(field-ball(clean)) + norm(field)/10**6,upper=True)
        true_distance = safe_float(norm(field-p)+rem,upper=True)
        if true_distance>=radius:
            raise ArithmeticError('The generating truth is outside the declared ball.')
        obs_id = f'compact_tm_g{int(gamma*100):03d}'
        np.savez_compressed(public/(obs_id+'.npz'),Y=y)
        records.append(dict(id=obs_id, object_id='compact_tm', gamma=gamma,
                            relative_noise=float(sample/np.linalg.norm(clean)),
                            sample_radius=sample, radius=radius, truth_distance_upper=true_distance))
    np.savez_compressed(public/'geometry.npz',G=geometry[0],S=geometry[1],B=geometry[2])
    np.savez_compressed(private/'reference.npz',material=truth,clean=clean,held=pristine_held,
                        held_sensor=held[1],held_incident=held[2])
    write(private/'generation.json',dict(design_landmark_distance=gap,
          numerical_landmark_port_error_upper=safe_float(legal_rem,upper=True),
          native_truth_peak=float(truth.max()),native_truth_mean=float(truth.mean()),
          inside_cells=int(np.count_nonzero(truth)),reference_train=audit,reference_held=held_audit,
          records=records,**elapsed(start)))
    write(public/'manifest.json',dict(physical_contract='SYNTHETIC_DESIGN.md',
          observations=[{k:v for k,v in row.items() if k not in ('truth_distance_upper',)} for row in records]))


def get_observations():
    public = OUT/'public'
    with np.load(public/'geometry.npz',allow_pickle=False) as saved:
        geometry = saved['G'], saved['S'], saved['B']
    rows = json.loads((public/'manifest.json').read_text())['observations']
    data = []
    for row in rows:
        with np.load(public/(row['id']+'.npz'),allow_pickle=False) as saved:
            data.append(Observation(row['id'],row['object_id'],100*row['relative_noise'],
                         row['sample_radius'],row['radius'],saved['Y']))
    return geometry, data


def decisions(geometry, observations):
    peak = np.linalg.norm(geometry[1],ord=2)**2
    results = []
    for obs in observations:
        model = Prepared(*geometry,obs.y)
        trials = []
        for scale in WEIGHT_SCALES:
            proof=certify(model, u=LOW, mu=scale*peak, radii={'fro':obs.padded_radius},
                          precision=112,threads=2)
            trials.append(dict(scale=scale, **proof))
        best = max(trials,key=lambda r:r.get('certified_frobenius_residual_lower',0.))
        d = best['certified_frobenius_residual_lower']
        result = dict(id=obs.id,cap=LOW,radius=obs.padded_radius,
                      source='fixed real-sector weight menu',
                      source_scale=best['scale'],strict_distance_lower=d,strict_ratio=d/obs.padded_radius,
                      strict_exclusion=d>obs.padded_radius,trials=trials,
                      requested_window_impossible=d>2.5*obs.padded_radius)
        results.append(result)
        print(obs.id, 'strict ratio', result['strict_ratio'], flush=True)
    return results


def assess(geometry, observations):
    f = OUT/'offline_truth'
    with np.load(f/'reference.npz',allow_pickle=False) as saved:
        truth, clean, held_clean = saved['material'], saved['clean'], saved['held']
        held = geometry[0], saved['held_sensor'],saved['held_incident']
    records = {}
    inside = truth>0
    for obs in observations:
        contrasts = {}
        name = OUT/'fits'/(obs.id+'.npz')
        receipt = json.loads(name.with_suffix('.json').read_text())
        with np.load(name,allow_pickle=False) as predictions:
            for policy in receipt['reports']:
                material = predictions[policy+'_material']
                p,rem,a = enclosure(geometry,material)
                held_p,held_r,audit_h = enclosure(held,material)
                r0 = norm(p-ball(obs.y))
                error = material-truth
                contrasts[policy] = dict(relative_l2=float(np.linalg.norm(error)/np.linalg.norm(truth)),
                    inside_relative_l2=float(np.linalg.norm(error[inside])/np.linalg.norm(truth[inside])),
                    background_rms=float(np.sqrt(np.mean(error[~inside]**2))),
                    train_clean_relative_l2=float(norm(p-ball(clean))/norm(ball(clean))),
                    withheld_relative_l2=float(norm(held_p-ball(held_clean))/norm(ball(held_clean))),
                    peak=float(material.max()), predicted_ports_recheck=a,predicted_held_recheck=audit_h,
                    train_residual_ratio_lower=max(0.,safe_float((r0-rem)/obs.padded_radius,upper=False)),
                    train_residual_ratio_upper=safe_float((r0+rem)/obs.padded_radius,upper=True))
        records[obs.id]=contrasts
    write(OUT/'offline_scores.json',records)



def linear_observability(geometry):
    _, sensor, incident = geometry
    kernel = np.einsum('ri,is->rsi', sensor, incident).reshape(-1, incident.shape[0])
    real_kernel = np.vstack((kernel.real, kernel.imag))
    sv = np.linalg.svd(real_kernel, compute_uv=False)
    fraction = sv**2/np.sum(sv**2)
    write(OUT/'offline_truth/linear_diagnostic.json', dict(
        scope='real contrast to stacked real and imaginary Born ports',
        singular_values=sv.tolist(), first_mode_squared_energy_fraction=float(fraction[0]),
        first_two_squared_energy_fraction=float(fraction[:2].sum())))


def main():
    ctx.prec,ctx.threads = 112,2
    observe()
    geometry,observations = get_observations()
    gates = decisions(geometry,observations)
    write(OUT/'gate.json',gates)
    native.SIDE,native.N,native.Model = NATIVE,NATIVE**2,DenseModel
    for obs,result in zip(observations,gates,strict=True):
        run(geometry,obs,result,OUT/'fits',stages=STAGES,released_cap=HIGH)
    assess(geometry,observations)
    linear_observability(geometry)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
