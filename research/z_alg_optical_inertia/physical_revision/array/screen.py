"""Frozen subspace screening. These scores never constitute native certificates."""

import json
from pathlib import Path
import sys
from time import process_time, perf_counter

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, svdvals
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
CURVE = STUDY / "submission/curve"
sys.path.insert(0, str(CURVE))
from search import MU, RATIOS, Projected, h, save_json, times

SOURCES = (10, 16, 24, 32)
RECEIVERS = (20, 32)
CAPS = (.75, 1., 1.25, 1.5, 2.)
CUTOFFS = (.1, .01, .001, 1.e-6)
GLOBAL_GRID = [(float(mu), float(t)) for t in RATIOS for mu in MU]


def load_json(path):
    return json.loads(path.read_text())


def frozen_profiles():
    source = load_json(CURVE / "global.75.json")["profiles"]
    quadrant = load_json(CURVE / "local.75.json")["profiles"][0]
    controls = [
        ("stable_scalar_optical", source[0]), ("scalar_real_01", source[1]),
        ("scalar_optical_05", source[5]), ("quadrant_optical_04", quadrant),
    ]
    for name, p in controls:
        yield dict(name=name, mu=p["scale"], partition=int(np.sqrt(len(p["weights"]))),
                   sector=p["weights"], optical=p["optical"], normalized=False)
    single = CURVE.parent / "single"
    profiles = load_json(single / "coherent_general_profiles.json")
    verified = load_json(single / "coherent_general_verified.json")["cases"]
    for name in sorted(verified):
        p = profiles[verified[name]["profile_index"]]
        yield dict(name=f"block4_{name[9:11]}", mu=p["mu"], partition=4,
                   sector=np.ravel(p["sector"]).tolist(), optical=np.ravel(p["optical"]).tolist(),
                   normalized=True)


def grid_profile(index):
    mu,t = GLOBAL_GRID[index]
    return dict(name=f"scalar_grid_{index}", mu=mu, partition=2,
                sector=[mu]*4, optical=[mu*t]*4, normalized=False)


def inputs():
    base = HERE / "observations"
    with np.load(STUDY / "synthetic/public/fine_geometry.npz", allow_pickle=False) as f:
        spec = f["spectrum"]
    with np.load(base / "arrays.npz", allow_pickle=False) as f:
        sensor = {p:f[f"S{p}"] for p in RECEIVERS}
        source = {s:f[f"B{s}"] for s in SOURCES}
    data = {p:{s:[] for s in SOURCES} for p in RECEIVERS}
    for x in load_json(base / "manifest.json")["observations"]:
        with np.load(base / x["file"], allow_pickle=False) as f:
            y = f["Y"]
        data[x["receivers"]][x["sources"]].append((x,y))
    return spec, sensor, source, data


def ranks(a):
    sig = svdvals(a, check_finite=False)
    relative = sig/sig[0]
    return dict(singular_values=sig.tolist(), relative_singular_values=relative.tolist(),
                effective_ranks={str(k):int(np.count_nonzero(relative>=k)) for k in CUTOFFS},
                numerical_rank=int(np.count_nonzero(relative>max(a.shape)*np.finfo(float).eps)),
                relative_min=float(relative[-1]), frobenius_norm=float(np.linalg.norm(a)))


def tails(hs, nu):
    spectra = np.linalg.eigvalsh((hs+hs.swapaxes(-1,-2).conj())/2)[...,::-1]
    return spectra, np.sqrt(np.maximum(spectra[...,nu:],0.).sum(axis=-1))


def frozen_run(proj, profile, u, groups):
    mu, partition = profile["mu"], profile["partition"]
    wr, opt = np.asarray(profile["sector"]), np.asarray(profile["optical"])
    if not profile["normalized"]:
        wr, opt = wr/mu, opt/mu
    d = proj["blocks"][partition]
    a = h(proj["S"].conj().T@proj["S"])/mu
    b = np.zeros_like(proj["B"])
    for wi, ti, (identity, real, imaginary, incident) in zip(wr,opt,d,strict=True):
        a += wi*(identity-u*real)+ti*imaginary
        b += (u*wi+1j*ti)/2*incident
    val, vec = eigh(h(a), check_finite=False, driver="evd")
    nu = int(np.count_nonzero(val<0.))
    responses = {}
    for s, gr in groups.items():
        if nu >= s or min(abs(val)) < 1.e-12:
            responses[s] = (np.zeros((len(gr["rows"]),s)),np.zeros(len(gr["rows"])))
            continue
        l = proj["S"].conj().T@gr["stack"]/mu
        l += np.tile(b[:,gr["idx"]],(1,len(gr["rows"])))
        v = (vec.conj().T@l).reshape(len(val),len(gr["rows"]),s)
        hs = gr["ygram"]-mu*np.einsum("kij,kil,k->ijl",v.conj(),v,1/val)
        responses[s] = tails(hs,nu)
    return nu, responses


def groups_for(data):
    groups={}
    left=0
    for s, entries in data.items():
        y = np.stack([j for _,j in entries])
        groups[s] = dict(rows=[r for r,_ in entries], y=y,
                         stack=np.concatenate(list(y),axis=1),
                         ygram=y.swapaxes(-2,-1).conj()@y,
                         idx=slice(left,left+s))
        left += s
    return groups


def projection(spec,sensor,incident):
    proj = Projected(spec,sensor,incident,7)
    q = proj.basis
    side = int(np.sqrt(q.shape[0]))
    t = fft2(q.reshape(side,side,-1),s=spec.shape,axes=(0,1))
    gq = ifft2(spec[:,:,None]*t,axes=(0,1))[:side,:side].reshape(q.shape)
    ij = np.arange(side)
    partition = 4
    labels = ((partition*ij[:,None]//side)*partition +
              partition*ij[None,:]//side).ravel()
    blocks=[]
    for group in range(partition**2):
        ind = labels == group
        u = q[ind].conj().T
        g = u@gq[ind]
        blocks.append((h(u@q[ind]),h(g),h(g/1j),u@incident[ind]))
    return dict(Gh=proj.Gh, Ga=proj.Ga, S=proj.S, B=proj.B, blocks={2:proj.blocks,4:blocks},
                dimension=proj.dimension, green_columns=proj.cost["projected_green_columns"]+q.shape[1])


def full_grid(proj,u,groups):
    d0 = proj["dimension"]
    results = {s: np.zeros((len(GLOBAL_GRID),len(gr["rows"]))) for s,gr in groups.items()}
    retained = {s: np.zeros((len(GLOBAL_GRID),len(gr["rows"]))) for s,gr in groups.items()}
    negative = np.full(len(GLOBAL_GRID),-1)
    for it,t in enumerate(RATIOS):
        ev,v = eigh(h(np.eye(d0)-u*proj["Gh"]+t*proj["Ga"]),check_finite=False,driver="evd")
        if min(abs(ev)) < 1.e-12:
            continue
        n0 = np.count_nonzero(ev<0.)
        receiver = proj["S"]@v
        dh = (v.conj().T@proj["B"])*(.5*(u+1j*t))
        sg = h((receiver/ev)@receiver.conj().T)
        g, w = eigh(sg,check_finite=False)
        q = (receiver/ev)@dh
        wd = w.conj().T
        for s,gr in groups.items():
            dd = dh[:,gr["idx"]]
            W = h((dd.conj().T/ev)@dd)
            deltas = wd@(gr["y"]-q[:,gr["idx"]])
            matrices = np.moveaxis(deltas.conj(),-2,-1)
            for im,mu in enumerate(MU):
                i = it*len(MU)+im
                gmu = g+mu
                if min(abs(gmu)) < mu*1.e-10:
                    continue
                nu = int(n0-np.count_nonzero(gmu<0.))
                if nu < 0:
                    continue
                negative[i] = nu
                if nu >= s:
                    continue
                hs = mu*((matrices/gmu)@deltas-W)
                modes, ds = tails(hs,nu)
                results[s][i] = ds
                retained[s][i] = modes[:,nu]
    return negative,results,retained


def write_rows(fd, profiles, groups, u, p, negative, result, retained, fixed):
    result_rows=[]
    for s,gr in groups.items():
        for column,meta in enumerate(gr["rows"]):
            score = result[s][:,column]
            order = sorted(np.flatnonzero((negative>=0)&(negative<s)),
                           key=lambda i:(score[i],retained[s][i,column]),reverse=True)[:8]
            row = dict(**meta,cap=u,min_projected_scalar_inertia=int(negative[negative>=0].min()),
                       scalar_grid_best=[dict(index=int(i),nu=int(negative[i]),distance=float(score[i]),
                                               first_retained=float(retained[s][i,column])) for i in order])
            row["fixed_profiles"]=[]
            for prof,(nu,rs) in zip(profiles,fixed,strict=True):
                e,v = rs[s]
                row["fixed_profiles"].append(dict(name=prof["name"],nu=nu,
                    distance=float(v[column]),first_retained=float(e[column,nu]) if nu<s else None,
                    source_eigenvalues=e[column].tolist() if nu<s else None))
            result_rows.append(row)
        fd[f"d_{p}_{u}_{s}"] = result[s]
        fd[f"eig_{p}_{u}_{s}"] = retained[s]
    fd[f"nu_{p}_{u}"] = negative
    return result_rows


def main():
    start, cpu = perf_counter(), process_time()
    spec,sensors,sources,observations = inputs()
    profiles=list(frozen_profiles())
    out=dict(scope="Numerical subspace proposals only", profiles=profiles,
             receivers={p:ranks(a) for p,a in sensors.items()},
             sources={s:ranks(a) for s,a in sources.items()},
             full_wave_number=2*np.pi*400e6/3e8, enclosing_radius_m=np.sqrt(2.),
             global_grid=GLOBAL_GRID, records=[], shared_projections={}, timings={})
    fd={}
    for p in RECEIVERS:
        t,c = perf_counter(),process_time()
        proj = projection(spec,sensors[p],np.column_stack(list(sources.values())))
        out["shared_projections"][p] = dict(**times(t,c), dimension=proj["dimension"],
                                             green_columns=proj["green_columns"])
        groups = groups_for(observations[p])
        for u in CAPS:
            t,c = perf_counter(),process_time()
            gs = full_grid(proj,u,groups)
            fs = [frozen_run(proj,q,u,groups) for q in profiles]
            rows = write_rows(fd,profiles,groups,u,p,*gs,fs)
            out["records"] += rows
            out["timings"][f"{p}/{u}"] = times(t,c)
            out.update(times(start,cpu))
            save_json(HERE/"screen.json",out)
            np.savez_compressed(HERE/"screen_scalar_grid.npz",**fd)
            print(p,u,[(s,sum(bool(x["scalar_grid_best"]) and x["scalar_grid_best"][0]["distance"]>x["padded_radius"]
                              for x in rows if x["sources"]==s),
                         sum(max(r["distance"] for r in x["fixed_profiles"])>x["padded_radius"]
                             for x in rows if x["sources"]==s))
                        for s in SOURCES],out["timings"][f"{p}/{u}"],flush=True)


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
