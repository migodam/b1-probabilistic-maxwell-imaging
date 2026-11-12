"""Generator-only new true3D cohort; online adapters never import this file."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,sys,time,argparse
from datetime import datetime,timezone
import numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.vector_adapter import DenseVectorModel
# The frozen R12 generator imports its sibling adapter by module name.
sys.path.insert(0,str(REPO/'research/research_max_struct_r12/delegated/vector_high'))
from research.research_max_struct_r12.delegated.vector_high.generate_high import volume_average


def dump(path,value):
    def cv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,default=cv,allow_nan=False)+'\n')


def geometries():
    def box(center,half,chi):return dict(center=np.array(center),half=np.array(half),chi=float(chi))
    # Independent analytic composites; last piece wins; coordinates in lambda.
    return dict(
      slotted_high_frame=[box([0,0,0],[.082,.075,.065],8),
          box([.006,.004,0],[.045,.038,.08],0),
          box([.055,.055,.025],[.038,.032,.018],0)],
      unequal_coupled_steps=[box([-.046,-.026,-.014],[.041,.047,.040],8),
          box([.049,.035,.025],[.036,.031,.053],3),
          box([.002,.007,.002],[.038,.017,.018],5)],
      mixed_crossed_bracket=[box([-.014,-.037,0],[.075,.029,.058],2),
          box([-.055,.012,.012],[.026,.074,.042],8),
          box([.034,.044,-.018],[.048,.021,.024],5),
          box([-.051,-.034,.020],[.018,.021,.023],0)])


def generate(out=HERE/'new3d'):
    out=Path(out)
    if out.exists():raise FileExistsError(f'preserve generated cohort: {out}')
    out.mkdir(parents=True);(out/'data_public').mkdir();(out/'sealed_truth').mkdir()
    names=list(geometries());cases=names+[names[0]+'_crossgrid8']
    contract=dict(status='FROZEN_BEFORE_OBSERVATION_GENERATION',created_utc=datetime.now(timezone.utc).isoformat(),
        cases=cases,independent_geometries=3,crossgrid_pair_not_independent=True,
        inverse_side=6,forward_side=[6,6,6,8],width_lambda=.24,material_box=[0,8],
        material='one shared real isotropic chi per cell, 3 vector current components',
        seed=20261113031,noise='iid CN(0,sigma²), sigma=.01||cleanY||_F/sqrt(number_complex_channels)',
        sources=60,receiver_vector_channels=180,acquisition='unchanged R8 three circles of electric dipoles; same physical positions across grids',
        initial_chi=.25,prior='unchanged R8 2e-4 edge mean square + 1e-6 material mean square',
        cell_phase_max={'inverse6':float(2*np.pi*.04*3),'forward8':float(2*np.pi*.03*3)},
        qualification='CM+RR finite point DDA, reduced internal cell phase; independent Mie control and crossgrid evidence required; no universal continuum certificate',
        truth_policy='generator-only truth is saved separately; online inputs contain Y,sigma,grid metadata only; all predictions must be sealed before offline scoring')
    dump(out/'FROZEN_PROTOCOL.json',contract)
    parts=geometries();dump(out/'sealed_truth/GEOMETRIES.json',parts)
    started,cpu=time.perf_counter(),time.process_time();coarse=DenseVectorModel();fine=None;rows=[];coarse_clean={}
    for i,cid in enumerate(cases):
        cross=cid.endswith('_crossgrid8');name=names[0] if cross else cid
        if cross and fine is None:fine=DenseVectorModel(side=8)
        model=fine if cross else coarse
        chi,integral=volume_average(model,parts[name],width=.24)
        inv_chi,inv_integral=volume_average(coarse,parts[name],width=.24)
        before=dict(model.cost);clean=model.forward(chi);st=model.state(chi)
        costs={k:v-before.get(k,0) for k,v in model.cost.items()}
        if not cross:coarse_clean[name]=clean.copy()
        rng=np.random.default_rng(contract['seed']+i)
        sigma=.01*np.linalg.norm(clean)/np.sqrt(clean.size)
        noise=sigma*(rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape))/np.sqrt(2)
        np.savez_compressed(out/'data_public'/f'{cid}.npz',Y=clean+noise,sigma=sigma,
            grid_shape=np.array([6,6,6]),width=np.array(.24),material_bounds=np.array([0,8]))
        np.savez_compressed(out/'sealed_truth'/f'{cid}.npz',chi=inv_chi,generation_chi=chi,
            generation_side=model.side,Y_clean=clean,noise=noise)
        rows.append(dict(case_id=cid,independent_geometry=name,generation_side=model.side,
            inverse_side=6,observation_shape=list(clean.shape),sigma=float(sigma),
            physical_costs=costs,forward_state_relative_residual=st['state_relative_residual'],
            generation_integral=integral,inverse_integral=inv_integral,
            fine_vs_coarse_clean_relative_difference=float(np.linalg.norm(clean-coarse_clean[name])/np.linalg.norm(clean)) if cross else None))
    dump(out/'GENERATION_RECEIPT.json',dict(status='COMPLETE',rows=rows,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started,
        inverse_setup=coarse.metadata,fine_setup=None if fine is None else fine.metadata))
    print(json.dumps(dict(status='GENERATED',directory=str(out),cases=cases)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=HERE/'new3d');a=p.parse_args();generate(a.out)
