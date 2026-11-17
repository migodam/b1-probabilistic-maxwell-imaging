"""+1% relative contrast on every nonactive host cell with original three physical shots."""
import json
import time
import numpy as np
from numpy.polynomial.polynomial import polyval
from patch_pencil import focused_polynomial, noisy_components
from pencil_experiment import CENTER, SCALE, ROUNDING, ROOT, REGIONS, host_model, ReducedKernel


def unpack(z):
    return np.asarray(z["real"])+1j*np.asarray(z["imag"])


def invert(d, model_R, tpl, q, radius):
    target=REGIONS["square"]
    g=focused_polynomial(d,model_R[target],tpl,CENTER,SCALE)
    components=[[CENTER+SCALE*l,CENTER+SCALE*h] for l,h in noisy_components(g,tpl["coord"],radius)]
    est=(components[0][0]+components[-1][-1])/2 if components else None
    t=(q-CENTER)/SCALE
    distance=abs(polyval(t,g))/np.linalg.norm(polyval(t,tpl["coord"].T))
    return dict(components=components,contains_truth=any(l<=q<=u for l,u in components),
                truth_profile_residual=float(distance),hull_midpoint=est,
                abs_midpoint_error=abs(est-q) if est is not None else None)


def metrics(rows):
    n=len(rows)
    finite=[r["abs_midpoint_error"] for r in rows if r["abs_midpoint_error"] is not None]
    return dict(n=n,covered=sum(r["contains_truth"] for r in rows),empty=sum(not r["components"] for r in rows),
                worst_midpoint_error_nonempty=max(finite,default=None),
                median_midpoint_error_nonempty=float(np.median(finite)) if finite else None,
                mismatch_l2_median=float(np.median([r["mismatch_l2"] for r in rows])),
                mismatch_l2_max=max(r["mismatch_l2"] for r in rows),
                mismatch_alone_exceeds_noise=sum(r["mismatch_l2"]>2e-4+ROUNDING for r in rows),
                mismatch_plus_noise_exceeds_noise=sum(r["total_l2"]>2e-4+ROUNDING for r in rows),
                truth_projection_exceeds_noise=sum(r["truth_profile_residual"]>2e-4+ROUNDING for r in rows))


def run():
    start=time.perf_counter()
    input_dir=ROOT/"results"
    raw=json.loads((input_dir/"pencil_raw.json").read_text())
    saved=np.load(input_dir/"pencil_observations.npz")
    out={}
    rows=[]
    vectors={}
    epsilon=2.e-4
    factor=1.01
    for layout in ["primary","confirmation"]:
        model,ker,_,_=host_model(layout)
        ref=ker.m0
        actual=factor*ref
        changed=ReducedKernel(model,actual,active=ker.active)
        empty_actual=model.forward(actual,independent=True)["y"].ravel()
        empty_ref=ker.y0.ravel()
        prefix=layout+"_square"
        design=raw["design"][layout]["regions"]["square"]["template"]
        tpl={"templates":unpack(design["templates"]),
             "coord":unpack(design["coord"])/np.sqrt(raw["common_total_energy"]/3.)}
        parameters=saved[prefix+"_truth_parameters"]
        rtarget=REGIONS["square"]
        nuisance=np.setdiff1d(np.arange(ker.r),rtarget)
        arms=["nominal","wrong_host_baseline","actual_host_baseline"]
        output={arm:[] for arm in arms}
        A=saved[prefix+"_00_pencil_feed"]
        out[layout]=dict(host_contrast_relative_l2=float(np.linalg.norm(actual-ref)/np.linalg.norm(ref)),
            max_abs_host_material_change=float(max(abs(actual-ref))),
            empty_host_three_shot_delta_norm=float(np.linalg.norm((empty_actual-empty_ref)@A)),
            H_rel=float(np.linalg.norm(changed.K-ker.K)/np.linalg.norm(ker.K)),
            B_rel=float(np.linalg.norm(changed.q-ker.q)/np.linalg.norm(ker.q)),
            R_rel=float(np.linalg.norm(changed.L-ker.L)/np.linalg.norm(ker.L)))
        for case,x in enumerate(parameters):
            name=prefix+f"_{case:02d}_pencil"
            A=saved[name+"_feed"]
            nominal_clean=saved[name+"_noise0"]
            nominal_noisy=saved[name+"_noise2"]
            noise=nominal_noisy-nominal_clean
            z=np.full(ker.r,x[0])
            z[nuisance]=x[1:]
            actual_y=model.forward(actual+ker.U@z,independent=True)["y"].ravel()
            clean=[nominal_clean,(actual_y-empty_ref)@A,(actual_y-empty_actual)@A]
            vectors[name+"_clean_rows"]=np.array(clean)
            vectors[name+"_absolute_measured"]=actual_y@A
            vectors[name+"_saved_noise"]=noise
            vectors[name+"_z"]=z
            for arm,y in zip(arms,clean):
                d=y+noise
                res=invert(d,ker.L[0],tpl,x[0],epsilon+ROUNDING)
                row=dict(layout=layout,case=case,arm=arm,truth=x[0],noise_l2=float(np.linalg.norm(noise)),
                         mismatch_l2=float(np.linalg.norm(y-nominal_clean)),
                         total_l2=float(np.linalg.norm(d-nominal_clean)),**res)
                output[arm].append(row)
                rows.append(row)
        out[layout]["comparison"]={arm:metrics(rs) for arm,rs in output.items()}
    result=dict(perturbation="Multiply all 24 nonactive contrasts by 1.01, with archive geometry, fixed square active material and archived three physical feeds.",
                radius=epsilon+ROUNDING,noise="Archived noise2 minus noise0 for each existing object. Exactly measured empty-host baseline is noiseless.",
                arms=arms,per_layout=out,combined={arm:metrics([r for r in rows if r["arm"]==arm]) for arm in arms},
                individual=rows,wall_seconds=time.perf_counter()-start)
    np.savez_compressed(input_dir/"geometry_mismatch_raw.npz",**vectors)
    (input_dir/"geometry_mismatch.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":run()
