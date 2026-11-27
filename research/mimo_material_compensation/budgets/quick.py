"""Compute preregistered uncompensated response curves."""
from pathlib import Path
import json
import os
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
sys.dont_write_bytecode = True

import numpy as np
from scipy.linalg import lu_factor, lu_solve

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1] / "z_alg_optical_inertia" / "blind_uq_mimo_theory"
CFG = json.loads((HERE / "freeze.json").read_text())
for folder in ("ports", "imaging", "diagnostics"):
    sys.path.insert(0, str(STUDY / folder))
from aperture_control import dark_apertures
from shared import Model, write_json, elapsed, now
from scipy.ndimage import label


def load_case(key):
    side = 96 if key.endswith("96") else 32
    case = f"n{side}_c02"
    with np.load(STUDY / "ports" / "public" / "operators" / f"n{side}_f400_a08.npz", allow_pickle=False) as op:
        base = {name: op[name] for name in op.files}
    _, exc, rec = dark_apertures(base)
    with np.load(STUDY / "results" / "extended" / "n32_c02.npz", allow_pickle=False) as obj:
        roi32 = obj["roi"].reshape(32, 32).astype(bool)
    if key=="joint32":
        with np.load(STUDY / "results" / "aperture_joint" / "n32_c02.npz", allow_pickle=False) as obj:
            center = obj["roi_zero_x"].ravel()
        roi = roi32.ravel()
    elif key=="native96":
        with np.load(STUDY / "results" / "gn_profiles" / "n96_c02" / "roi_02.npz", allow_pickle=False) as obj:
            center = obj["material"].ravel()
            roi = obj["roi"].ravel().astype(bool)
    elif key=="air32":
        center = np.zeros(side**2)
        roi = roi32.ravel()
    elif key=="air96":
        center = np.zeros(side**2)
        roi = roi32.repeat(3, axis=0).repeat(3, axis=1).ravel()
    if np.any(center[roi] != 0) or np.any(center<0) or np.any(center>4.5):
        raise ValueError("Invalid public reference and query mask")
    observed = {}
    with np.load(STUDY / "ports" / "public" / "cases" / f"{case}.npz", allow_pickle=False) as obj:
        radius = float(obj["radius"])
        observation_id = str(obj["observation_id"].item())
        if key in ("joint32", "native96"):
            observed["A"] = obj["Y"]
    if key=="joint32":
        with np.load(STUDY / "results" / "aperture_control" / f"{case}.npz", allow_pickle=False) as obj:
            observed["B"] = obj["B_interleaved_Y"]
    return side, base, exc, rec, roi, center, radius, observation_id, observed


class Response:
    def __init__(self, side, base, exc, rec):
        self.base, self.exc, self.rec = base, exc, rec
        self.g = None
        self.model = Model(dict(base, S=rec, E0=exc), np.zeros((len(rec),exc.shape[1]), complex),
                           1., "fixed-ports", side)
        self.last = {}

    def full(self, x, derivatives=False):
        m=self.model
        if np.array_equal(x, self.last.get("x")) and ("u" in self.last or not derivatives):
            return self.last
        ids = np.flatnonzero(x)
        if len(ids)<=1000:
            # All inactive currents vanish exactly. Numerical proposals can use
            # the independently checked released Fourier Green blocks.
            from scipy.fft import ifft2
            side = m.side
            yr,xr = divmod(ids,side)
            lag = ifft2(self.base["spectrum"])
            d = lag.shape[0]
            g = lag[(yr[:,None]-yr[None,:])%d,(xr[:,None]-xr[None,:])%d]
            g += self.base["low_rank_left"][ids] @ self.base["low_rank_right"][:,ids]
            j = np.zeros_like(self.exc)
            if len(ids):
                factor = lu_factor(np.eye(len(ids))-x[ids,None]*g, check_finite=False)
                j[ids] = lu_solve(factor, x[ids,None]*self.exc[ids], check_finite=False)
                m.cost["physical_lu"] += 1
            m.cost["physical_forward_rhs"] += self.exc.shape[1]
            if derivatives:
                if len(ids):
                    # W A=C implies W_I=C_I A_II^-1 and W=C+W_I X_I G_I,:.
                    adjoint_current = np.zeros((len(x),len(self.rec)), complex)
                    u_active = lu_solve(factor,self.rec[:,ids].T,trans=1,check_finite=False).T
                    adjoint_current[ids]=(u_active*x[ids]).conj().T
                    u = self.rec + m.green.adjoint(adjoint_current).conj().T
                    m.cost["physical_adjoint_rhs"] += len(self.rec)
                else:
                    u = self.rec.copy()
        else:
            j = m.solve(x,x[:,None]*self.exc)
            if derivatives:
                u = m.solve(x,self.rec.conj().T,adjoint=True).conj().T
        p=self.rec@j
        e=self.exc+m.green.apply(j)
        self.last = dict(x=x.copy(),j=j,e=e,p=p)
        if derivatives:
            self.last["u"]=u
        return self.last


def block(a, key):
    layout = CFG["layout"][key]
    return a[np.ix_(layout["rx"],layout["tx"])]


def main(keys):
    for key in keys:
        start=now()
        side, base, exc, rec, roi, center, r, obs_id, observed=load_case(key)
        response=Response(side,base,exc,rec)
        first=response.full(center)["p"]
        records=[]
        for t in CFG["amplitudes"]:
            second=response.full(center+t*roi)["p"]
            results={}
            for layout in ("A","B"):
                a=block(first,layout)
                b=block(second,layout)
                results[layout]={"prospective":float(np.linalg.norm(b-a)/(2*r))}
                if layout in observed:
                    results[layout].update(center_over_r=float(np.linalg.norm(a-observed[layout])/r),
                                           actual_over_r=float(np.linalg.norm(b-observed[layout])/r))
            records.append(dict(t=t,**results))
        doc=dict(key=key, side=side, cells=int(roi.sum()), areas_m2=float(roi.sum())*(2/side)**2,
                 components=label(roi.reshape(side,side))[1], r=r, observation_id=obs_id,
                 **elapsed(start), cost=dict(response.model.cost), records=records)
        write_json(HERE/f"quick_{key}.json",doc)
        print(key, [(i["t"],round(i["A"]["prospective"],2),round(i["B"]["prospective"],2),
                     None if "actual_over_r" not in i["A"] else round(i["A"]["actual_over_r"],2),
                     None if "actual_over_r" not in i["B"] else round(i["B"]["actual_over_r"],2)) for i in records], flush=True)


if __name__=="__main__":
    main(sys.argv[1:] or list(CFG["reference_tasks"]))
