"""Global cancellation and available Fourier-current directions at early CIE iterates."""
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
SRAP = RESEARCH / 'srap_cie_array_study_20261121'
sys.path.insert(0, str(SRAP))
import engine
import run_study
from check_currents import imag_h, local_readout


def cases():
    entries = {c['case']: c for c in run_study.cases()}
    for name in ('la180_austria_rot90_eps3p5', 'optfresh_01__pw180_10x20__n02',
                 's2d_004__pw180_10x20', 'fresnel_twin_8GHz_half180_tx'):
        if name.startswith('optfresh'):
            folder = RESEARCH / 'z_alg_optical_inertia/synthetic/public'
            c = dict(entries['s2d_004__pw180_10x20'], case=name)
            p = run_study.load_public(entries['s2d_004__pw180_10x20'])
            geom = np.load(folder / 'pw180_10x20.npz')
            if not np.array_equal(p['S'], geom['S']) or not np.array_equal(p['E'], geom['B']):
                raise ValueError('Native archive and inversion array do not match')
            p['Y'] = np.load(folder / (name+'.npz'))['Y']
            contract = json.loads((folder/'manifest.json').read_text())
            item = next(v for v in contract['observations'] if v['id'] == name)
            p['delta'] = item['radius_noise_fro']
        else:
            c = entries[name]
            p = run_study.load_public(c)
        yield name, c, p


def similarity(x, y):
    return float(np.vdot(x,y).real / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-300))


def measure(model, y, m, j, e, u, beta, step):
    b = model.E
    vals = local_readout(j, e, b, model.side)
    gains = np.linalg.norm(b,axis=0) / np.sqrt(len(b))
    jw, ew, bw = j / gains, e / gains, b / gains
    q = imag_h(jw.conj().T @ ew)
    rj = ((ew-bw) - model.green.adjoint(jw)) / (2j)
    mult = rj - .5j*bw
    d = model.ports.receive(j)-y
    if model.ports.valid_mask is not None:
        d = np.where(model.ports.valid_mask, d, 0.)
    wd = 1 / np.maximum(np.sum(abs(y)**2, axis=0),1e-300)
    ws = 1 / np.sum(abs(b)**2, axis=0)
    p = (beta/(1+beta*m))[:,None]
    defect = p*p*(j-m[:,None]*e)*ws
    raw = model.ports.adjoint(d*wd) + defect - model.green.adjoint(m[:,None]*defect)
    baseline = u.adjoint(raw)
    for name, q_part in [('diagonal', np.diag(q.diagonal().real)), ('full',q)]:
        g = (mult @ q_part) / gains
        g_proj = u.adjoint(g)
        vals[name+'_gradient_projection_ratio'] = float(np.linalg.norm(g_proj)/max(np.linalg.norm(g),1e-300))
        vals[name+'_gradient_cosine_with_cie'] = similarity(g_proj, baseline)
    vals.update(step=step, n=len(b), s=b.shape[1], auxiliary_port_relative=float(np.linalg.norm(d)/np.linalg.norm(y)))
    return vals


def run():
    result = {}
    for name,c,p in cases():
        tic, cpu = perf_counter(), process_time()
        model = engine.legacy.PaperModel(p,c,engine.legacy.Clock(400))
        y = p['Y']
        iterator = engine.run_cie_arm(model,y,c,'dual',np.zeros(model.n),float(p['delta']))
        rows = []
        for m,step,stage,_ in iterator:
            if step in (1,16,64):
                a = iterator.gi_frame.f_locals
                if step==1:
                    rows.append(measure(model,y,np.zeros_like(m),a['JD'],a['QD'],a['U'],a['bj'],0))
                rows.append(measure(model,y,m,a['J'],a['Q'],a['U'],a['bj'],step))
            if step>=64:
                break
        iterator.close()
        result[name] = dict(measurements=rows, wall_s=perf_counter()-tic,cpu_s=process_time()-cpu,
                            ops=dict(model.clock.counts))
        (HERE/'early_directions.json').write_text(json.dumps(result,indent=2)+'\n')
        print(name,[(r['step'],round(r['coherent_fraction'],3),round(r['cost_of_sharing_fraction'],3),
                     round(r['full_gradient_cosine_with_cie'],3)) for r in rows],flush=True)


if __name__ == '__main__':
    with threadpool_limits(1): run()
