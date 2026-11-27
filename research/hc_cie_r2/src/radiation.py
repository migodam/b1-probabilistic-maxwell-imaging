"""Quantify the constant real weight using an independently sampled radiation frame."""
import argparse
from pathlib import Path
from time import perf_counter,process_time
from scipy.special import jv
import json
from spectra_inputs import load_observations
from interaction import LAB, np, NORM, threadpool_limits
from scipy.linalg import null_space,lu_factor,lu_solve


def evaluate(workspace):
    design=json.loads((LAB/"configs/spectra.json").read_text())
    design.update(r13_cases=["r13_fresh_0000","r13_fresh_0001"],r15_cases=[],r16_cases={},
                  srap_objects=["s2d_010"],srap_layouts=["pw360","ls_line_same"])
    rows=[]
    for obs in load_observations(workspace,design):
        wall,cpu=perf_counter(),process_time()
        with np.load(workspace/obs.provenance["observation"]) as data:
            x=data["coords"] if "coords" in data else data["centers"]
            k=float(data["k"]);rad=float(data["radius"])
        b=obs.b;g=obs.g;s=obs.s;m=obs.base.reshape(-1)
        h=(g-g.conj().T)/(2j)
        dark=null_space(s)
        jc=lu_solve(lu_factor(np.eye(m.size)-m[:,None]*g),m[:,None]*b)
        meas=NORM(dark.conj().T@(b+2j*h@jc),2)
        factor=float(.5*np.pi*k*rad*jv(1,k*rad))
        for n in [128,256]:
            angle=2*np.pi*np.arange(n)/n
            directions=np.stack([np.cos(angle),np.sin(angle)],axis=-1)
            c=np.sqrt(factor/n)*np.exp(-1j*k*(directions@x.T))
            leak=NORM(c@dark,2)
            err=NORM(h-c.conj().T@c,"fro")
            l,ss,r=np.linalg.svd(c.conj().T,full_matrices=False)
            keep=ss>leak
            aa=(r[keep].conj().T/ss[keep])@(l[:,keep].conj().T@b)
            illumination_def=NORM(b-c.conj().T@aa,2)
            bound=illumination_def+leak*(NORM(aa,2)+2*NORM(c@jc,2))+2*err*NORM(jc,2)
            rows.append(dict(observation=obs.key,directions=n,frame_error_fro=err,frame_scale=factor,
                       dark_radiation_operator=leak,represented_source_rank=int(keep.sum()),
                       source_coefficient_spectral_norm=NORM(aa,2),source_approximation_spectral_norm=illumination_def,
                       B_norm=NORM(b,2),measured_linear_norm=meas,linear_upper_physical=bound,
                       ratio_measured_to_B=meas/NORM(b,2),ratio_upper_to_B=bound/NORM(b,2),
                       unit_weight_scale=1/np.sqrt(m.size)))
        rows[-1].update(cpu_s=process_time()-cpu,wall_s=perf_counter()-wall,full_lu=1,forward_rhs=b.shape[1])
    return rows


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("workspace",type=Path);args=p.parse_args()
    with threadpool_limits(1):r=evaluate(args.workspace.resolve())
    (LAB/"results/real_weights/radiation.json").write_text(json.dumps(r,indent=2)+"\n")
    for row in r[1::2]: print(row["observation"],row["dark_radiation_operator"],row["ratio_measured_to_B"],row["ratio_upper_to_B"])
