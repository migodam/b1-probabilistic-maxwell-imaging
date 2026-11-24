"""Feasible witnesses for the complex-symmetric bounded transfer relaxation."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import svd
from flint import arb, ctx
from threadpoolctl import threadpool_limits

from branches import as_branch, expand
from contraction_control import spectral_upper, norm
from native_box import HERE, ball, clock, elapsed, hi
from site_driver import DEFAULT_PUBLIC


class ReciprocalControl:
    def __init__(self, sensor, incident, gain):
        self.sensor = sensor
        self.incident = incident
        self.gain = float(gain)
        q, scales, _ = svd(np.column_stack((sensor.conj().T, incident.conj())),
                           full_matrices=False, check_finite=False)
        keep = scales > scales[0]*(128*np.finfo(float).eps)
        self.q = q[:, keep]
        self.rank = q[:, keep].shape[1]
        i, j = self.i, self.j = np.triu_indices(self.rank)
        self.off = i != j
        self.scale = np.where(self.off, np.sqrt(2.), 1.)
        a, b = sensor @ self.q, self.q.T @ incident
        values = a[:, i][:, None, :]*b[j, :].T[None, :, :]
        values[..., self.off] += a[:, j[self.off]][:, None, :]*b[i[self.off], :].T[None, :, :]
        self.map = float(gain)*values.reshape(sensor.shape[0]*incident.shape[1], -1)/self.scale
        u, self.d, self.vh = svd(self.map, full_matrices=False, check_finite=False)
        self.u = u
        self.arb_forms = None

    def symmetric(self, x):
        matrix = np.empty((self.rank, self.rank), dtype=complex)
        upper = x/self.scale
        matrix[self.i, self.j] = upper
        matrix[self.j, self.i] = upper
        return matrix

    def coefficients(self, matrix):
        return self.scale*matrix[self.i, self.j]

    def spectral_projection(self, values):
        matrix = self.symmetric(values)
        u, d, vh = svd(matrix, full_matrices=False, check_finite=False)
        trial = (u*np.minimum(d, .99995)) @ vh
        return self.coefficients((trial+trial.T)*.5)

    def project_data(self, x, transformed_y, perpendicular, radius):
        # The minimum-change projection in the economy singular basis.
        e = transformed_y-self.d*(self.vh @ x)
        remaining = radius*radius-perpendicular*perpendicular
        if remaining <= 0:
            return None
        if np.vdot(e, e).real <= remaining:
            return x
        scale = self.d[0]
        relative_spectrum = self.d/scale
        squared = relative_spectrum*relative_spectrum

        def root(trial):
            component = e/(1+trial*squared)
            return np.vdot(component, component).real-remaining

        lower, upper = 0., 1.
        while root(upper) > 0 and upper < 1.e30:
            upper *= 16.
        if root(upper) > 0:
            return None
        for _ in range(48):
            middle = (lower+upper)/2.
            if root(middle) > 0:
                lower = middle
            else:
                upper = middle
        correction = upper*relative_spectrum*e/(scale*(1+upper*squared))
        return x+self.vh.conj().T @ correction

    def propose(self, y, radius, iterations):
        wave = y.ravel()
        transformed = self.u.conj().T @ wave
        perpendicular = float(np.linalg.norm(wave-self.u @ transformed))
        goal = .9992*radius
        x0 = np.zeros(self.map.shape[1], complex)
        z = self.project_data(x0, transformed, perpendicular, goal)
        maximum = float(np.linalg.norm(wave))
        best = x0
        if z is None:
            return self.symmetric(best), maximum, 0
        for step in range(1,iterations+1):
            x = self.spectral_projection(z)
            distance = float(np.linalg.norm(self.map @ x-wave))
            if distance<maximum:
                maximum, best = distance,x
            if maximum<.99998*radius:
                break
            reflected = self.project_data(2*x-z,transformed,perpendicular,goal)
            if reflected is None:
                break
            z += reflected-x
        return self.symmetric(best),maximum,step

    def verify(self, x, y, radius):
        bits = x.view(np.uint64).reshape(len(x),len(x),2)
        if not np.array_equal(bits,bits.transpose(1,0,2)):
            raise ArithmeticError("The stored matrix entries are not exactly symmetric.")
        if self.arb_forms is None:
            Q = ball(self.q)
            self.arb_forms = (ball(self.sensor)*Q, Q.transpose()*ball(self.incident), spectral_upper(Q)**2)
        left,right,q2 = self.arb_forms
        X = ball(x)
        bound = (q2*spectral_upper(X)).upper()
        defect = norm(ball(y)-arb(self.gain)*left*X*right).upper()
        return dict(norm_Q_squared_upper=hi(q2),norm_T_over_C_upper=hi(bound),
                    residual_upper=hi(defect),relative_radius_upper=hi(defect/arb(radius)),
                    matrix_entry_bits_symmetric=True,
                    strict_fit=bool(bound <= 1 and defect <= arb(radius)))


def run(public, contractions, iterations, output, witnesses):
    ctx.prec,ctx.threads=128,1
    started=clock()
    rows=json.loads(contractions.read_text())["rows"]
    cases={item["case"]:item for item in json.loads((public/"manifest.json").read_text())["cases"]}
    states={}
    geometry={}
    with np.load(public/"geometry96.npz",allow_pickle=False) as raw:
        for ap in sorted({c["array"] for c in cases.values()}):
            geometry[ap]={"S":raw[ap+"__S"],"B":raw[ap+"__B"]}
    data={}
    records,store=[],{}
    for old in rows:
        pair=clock()
        name, coordinates = old["observation"],tuple(old["tile"])
        case=cases[name]
        ap=case["array"]
        if name not in data:
            with np.load(public/case["measurement_file"],allow_pickle=False) as raw:
                data[name]=(raw["Y"],float(raw["delta"]))
        Y,r=data[name]
        if r!=old["radius"]:
            raise ValueError("The exposed comparison changed its noise radius.")
        signature=ap,coordinates,old["current_gain_C"]
        reused=signature in states
        if not reused:
            cap=np.zeros((4,4))
            cap[coordinates]=4.
            support=expand(as_branch("one_tile",cap),96)>0
            exact=geometry[ap]
            model=ReciprocalControl(exact["S"][:,support],exact["B"][support],signature[-1])
            q_key=f"{ap}__tile{coordinates[0]}_{coordinates[1]}__q{len(states)}"
            states[signature]=(model,q_key)
        model,q_key=states[signature]
        setup=elapsed(pair)
        then=clock()
        X,proposed,count=model.propose(Y,r,iterations)
        numeric=elapsed(then)
        row=dict(observation=name,tile=coordinates,original_decision=old["original_decision"],
                 arbitrary_contraction_fit=old["strict_energy_only_fit"],radius=r,current_gain_C=model.gain,
                 basis_rank=model.rank,observation_map=[int(a) for a in model.map.shape],
                 reused_setup=reused,iterations=count,numeric_residual=proposed,
                 setup=setup,proposal=numeric,strict_fit=False)
        if proposed<r:
            then=clock()
            row.update(model.verify(X,Y,r))
            row["proof"]=elapsed(then)
            if row["strict_fit"]:
                key=f"{name}__tile{coordinates[0]}_{coordinates[1]}__X"
                store[key]=X
                store[q_key]=model.q
                row.update(witness_matrix=key,witness_basis=q_key)
        row.update(elapsed(pair))
        records.append(row)
    with witnesses.open("wb") as fd:
        np.savez_compressed(fd,**store)
    wrong=[row for row in records if row["original_decision"]=="PRUNE"]
    results=dict(
        scope="Exposed compact-site observations; pure bounded reciprocal transfer. No material law.",
        method="Frobenius-orthonormal symmetric basis, economy SVD and at most the stated number of DR steps; Arb witnesses",
        max_iterations=iterations, precision_bits=ctx.prec,
        geometry_factors=len(states),input_pairs=len(records),verified_fits=sum(r["strict_fit"] for r in records),
        active_excluded_pairs=len(wrong),fits_in_active_excluded=sum(r["strict_fit"] for r in wrong),
        prior_unrestricted_fits_in_active_excluded=sum(r["arbitrary_contraction_fit"] for r in wrong),
        rows=records,**elapsed(started))
    output.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")
    print(len(records),results["verified_fits"],len(wrong),results["fits_in_active_excluded"],
          results["wall_seconds"],results["cpu_seconds"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public",type=Path,default=DEFAULT_PUBLIC)
    parser.add_argument("--contractions",type=Path,default=HERE/"contraction_control.json")
    parser.add_argument("--iterations",type=int,default=100)
    parser.add_argument("--output",type=Path,default=HERE/"reciprocal_control.json")
    parser.add_argument("--witnesses",type=Path,default=HERE/"reciprocal_witnesses.npz")
    args=parser.parse_args()
    with threadpool_limits(1):
        run(args.public,args.contractions,args.iterations,args.output,args.witnesses)
