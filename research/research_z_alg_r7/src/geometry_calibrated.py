"""One unloaded-fixture exposure for calibration of one scalar host contrast scale."""
import json
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from patch_pencil import design_pencil, whitening
from pencil_experiment import CENTER, SCALE, ROOT, REGIONS, ROUNDING, ReducedKernel, host_model, target_reciprocal
from geometry_mismatch import metrics, invert, unpack


EPS = 2e-4
CAL_SEED = 1811202603
SCALE_BOX = (.97,1.03)
SIMULATED_SCALE = 1.01
TARGET_ENERGY = 2.
FEED_CAP = 4.


def physical_read(model, material):
    return model.forward(material,independent=True)["y"].ravel()


def unit_scale_feed(model, nominal_host, gram):
    state=model.forward(nominal_host,independent=True)
    receive=np.linalg.solve(state["L"].T,model.S.T).T
    slope=((receive*nominal_host)@state["E"]).ravel()
    identity=np.eye(model.ns)
    def weighted(lam):
        a=np.linalg.solve(gram+lam*identity,slope.conj())
        return a/np.sqrt(np.vdot(a,gram@a).real)
    a=weighted(0.)
    if np.linalg.norm(a)>FEED_CAP:
        b=slope.conj()*FEED_CAP/np.linalg.norm(slope)
        if np.vdot(b,gram@b).real<=1:
            a=b
        else:
            lam=np.exp(brentq(lambda loglambda:np.linalg.norm(weighted(np.exp(loglambda)))-FEED_CAP,-35,35))
            a=weighted(lam)
    return a,slope@a


def fit_scale(model, nominal_host, read, beam):
    def objective(scale):
        return abs(physical_read(model,scale*nominal_host)@beam-read)**2
    estimate=minimize_scalar(objective,method="bounded",bounds=SCALE_BOX,options={"xatol":1e-14})
    return estimate.x, np.sqrt(estimate.fun), estimate.nfev


def redesign(model, host, active, W, saved_basis):
    ker=ReducedKernel(model,host,active=active)
    H,B,R=ker.K,ker.q,ker.L.ravel()
    vnom=target_reciprocal(H,R,np.full(ker.r,CENTER))
    target=REGIONS["square"]
    design=design_pencil(H,B@W,target,CENTER,SCALE,vnom)
    if design is None or design["shots"].shape[1]!=3:
        raise RuntimeError("Recalibrated square does not have a three-shot polynomial pencil")
    # Choose the orthonormal basis closest to the archived physical shots in the radiation metric.
    # This preserves the pencil and allows the exact same three-channel noise vectors to be used.
    P=design["shots"]
    u,_,vh=np.linalg.svd(P.conj().T@saved_basis)
    rotation=u@vh
    amplitude=np.sqrt(TARGET_ENERGY/3)
    feeds=amplitude*W@P@rotation
    coordinates=rotation.conj().T@design["coord"]/amplitude
    return ker,feeds,dict(templates=design["templates"],coord=coordinates),design["relative_synthesis"]


def run():
    start=time.perf_counter()
    folder=ROOT/"results"
    base=json.loads((folder/"pencil_raw.json").read_text())
    archives=np.load(folder/"pencil_observations.npz")
    rng=np.random.default_rng(CAL_SEED)
    raw,summary,rows={},{},[]
    for layout in ["primary","confirmation"]:
        model,nominal,gram,W=host_model(layout)
        actual_host=SIMULATED_SCALE*nominal.m0
        a,slope=unit_scale_feed(model,nominal.m0,gram)
        phase=rng.uniform(0,2*np.pi)
        cal_noise=EPS*np.exp(1j*phase)
        ycal=physical_read(model,actual_host)@a+cal_noise
        estimate,fit_res,nfev=fit_scale(model,nominal.m0,ycal,a)
        old=base["design"][layout]["regions"]["square"]["template"]
        basis=unpack(old["shots"])
        raw.update({layout+"_cal_feed":a,layout+"_cal_noise":cal_noise,layout+"_cal_reading":ycal,
                    layout+"_nominal_derivative":slope})
        calibration=dict(fitted_scale=estimate, scale_error=estimate-SIMULATED_SCALE,
                         cal_power=np.vdot(a,gram@a).real,cal_feed_norm=np.linalg.norm(a),
                         noise_phase_rad=phase,slope_abs=abs(slope),cal_residual=fit_res,objective_evaluations=nfev)
        layout_arms={}
        for mode,used_scale in [("one_calibration",estimate),("true_scale_diagnostic",SIMULATED_SCALE)]:
            used_host=used_scale*nominal.m0
            ker,feeds,tpl,synth=redesign(model,used_host,nominal.active,W,basis)
            host_pred=ker.y0.ravel()@feeds
            prefix=layout+"_square"
            xs=archives[prefix+"_truth_parameters"]
            nonpatch=np.setdiff1d(np.arange(nominal.r),REGIONS["square"])
            individual=[]
            for case,x in enumerate(xs):
                name=prefix+f"_{case:02d}_pencil"
                z=np.full(nominal.r,x[0])
                z[nonpatch]=x[1:]
                actual=physical_read(model,actual_host+nominal.U@z)@feeds
                clean=actual-host_pred
                # Used only for quantifying forward error after the fit, never for choosing beams/scale.
                matched_model=physical_read(model,used_host+nominal.U@z)@feeds-host_pred
                noise=archives[name+"_noise2"]-archives[name+"_noise0"]
                data=clean+noise
                result=invert(data,ker.L[0],tpl,x[0],EPS+ROUNDING)
                row=dict(layout=layout,mode=mode,case=case,truth=x[0],noise_l2=float(np.linalg.norm(noise)),
                         mismatch_l2=float(np.linalg.norm(clean-matched_model)),
                         total_l2=float(np.linalg.norm(data-matched_model)),**result)
                individual.append(row)
                rows.append(row)
                raw[name+"_"+mode+"_response"]=data
                raw[name+"_"+mode+"_noiseless_experiment"]=actual
                raw[name+"_"+mode+"_model_error"]=clean-matched_model
            total_energy=np.trace(feeds.conj().T@gram@feeds).real
            max_norm=max(np.linalg.norm(feeds,axis=0))
            if total_energy>TARGET_ENERGY+1e-10 or max_norm>FEED_CAP+1e-10:
                raise ValueError("Recalibrated physical shots exceeded the design power or feed caps")
            raw[layout+"_"+mode+"_feeds"]=feeds
            raw[layout+"_"+mode+"_predicted_empty"]=host_pred
            raw[layout+"_"+mode+"_coords"]=tpl["coord"]
            raw[layout+"_"+mode+"_templates"]=tpl["templates"]
            layout_arms[mode]=dict(**metrics(individual),used_scale=used_scale,physical_total_energy=total_energy,
                                   max_feed_norm=max_norm,relative_synthesis=synth)
        summary[layout]=dict(calibration=calibration,arms=layout_arms)
    result=dict(simulated_host_scale=SIMULATED_SCALE,allowed_scale_range=list(SCALE_BOX),
                common_complex_noise_radius=EPS,calibration_noise_seed=CAL_SEED,
                model="One unknown real global nonactive contrast scale; 1 unloaded-fixture shot reused across 16 objects for each layout.",
                error_rule="A single uniformly random calibration noise phase per layout with fixed magnitude 2e-4. Archived original complex noise vectors on redesigned three shots.",
                target_total_energy=TARGET_ENERGY,one_calibration_energy_max=1.,layouts=summary,
                combined={mode:metrics([r for r in rows if r["mode"]==mode]) for mode in
                          ["one_calibration","true_scale_diagnostic"]},individual=rows,wall_seconds=time.perf_counter()-start)
    (folder/"geometry_calibrated.json").write_text(json.dumps(result,indent=2)+"\n")
    np.savez_compressed(folder/"geometry_calibrated_raw.npz",**raw)


if __name__=="__main__":run()
