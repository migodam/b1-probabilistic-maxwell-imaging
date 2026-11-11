"""Freeze and generate two independent high-chi true3D geometries plus pair.

Material labels are available only to this generator and later offline scorer.
Volumetric truth uses exact integrals of piecewise-constant rectilinear solids.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from datetime import datetime,timezone
import sys,json,time
import numpy as np

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3]
sys.path.insert(0,str(REPO))
from model_adapter import HighContrastVectorPR3D


CASES=['random_high8','random_mixed2_8','random_high8_crossgrid8']
METHODS=['CIE_sourcewise','Dual_sourcewise','DR','AP']


def dump(path,value):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,default=convert,allow_nan=False)+'\n')


def geometries(seed):
    """Independent random draws, fixed before any observation or inversion."""
    a,b=[np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(2)]
    def box(rng,center,half,value):
        return dict(center=np.asarray(center)+rng.uniform(-.025,.025,3),
                    half=np.asarray(half)*rng.uniform(.86,1.14,3),chi=float(value))
    high=[box(a,[-.025,.005,-.02],[.16,.13,.15],8),
          box(a,[-.14,.085,.075],[.075,.07,.12],8),
          box(a,[.12,-.105,-.055],[.12,.065,.09],8),
          box(a,[.025,.025,.015],[.045,.04,.055],0)]
    mixed=[box(b,[.02,-.015,.005],[.18,.15,.16],2),
           box(b,[-.095,.06,-.055],[.09,.09,.105],8),
           box(b,[.11,-.095,.08],[.08,.065,.09],8),
           box(b,[-.005,.01,.045],[.045,.04,.055],0)]
    return {'random_high8':high,'random_mixed2_8':mixed}


def volume_average(model,parts,width=.6):
    """Exact integrate last-piece-wins boxes over every Maxwell grid cell.

    Splitting at all object planes creates disjoint constant-value bricks.
    Intersections with each inverse/forward cell are analytic rectangular
    volumes. Midpoints merely identify each constant brick's label.
    """
    planes=[]
    for axis in range(3):
        cuts=[-width/2,width/2]
        for part in parts:
            cuts.extend([float(part['center'][axis]-part['half'][axis]),
                         float(part['center'][axis]+part['half'][axis])])
        planes.append(np.unique(np.clip(cuts,-width/2,width/2)))
    answer=np.zeros(model.N);volume=0.;bricks=0
    for ix in np.ndindex(*(len(p)-1 for p in planes)):
        lo=np.array([p[i] for p,i in zip(planes,ix)])
        hi=np.array([p[i+1] for p,i in zip(planes,ix)])
        center=(lo+hi)/2;value=0.
        for part in parts:
            if np.all(abs(center-part['center'])<part['half']):value=part['chi']
        if value==0:continue
        overlap=np.maximum(np.minimum(model.points+model.h/2,hi)-np.maximum(model.points-model.h/2,lo),0)
        answer+=value*np.prod(overlap,axis=1)/model.h**3
        volume+=np.prod(hi-lo);bricks+=1
    answer=np.clip(answer,0,8)
    return answer,dict(nonzero_disjoint_bricks=bricks,support_volume=volume,
                       integrated_contrast=float(np.sum(answer)*model.h**3),
                       rule='analytic volume of disjoint constant-value bricks; no cell-center label approximation')


def main():
    if (HERE/'FROZEN_PROTOCOL.json').exists():raise FileExistsError('preserve frozen generation')
    protocol=dict(status='FROZEN_BEFORE_RANDOM_GEOMETRY_AND_NEW_OUTCOMES',
        created_utc=datetime.now(timezone.utc).isoformat(),seed=2026111213,
        cases=CASES,independent_geometries=2,crossgrid_pair_not_independent=True,methods=METHODS,
        inverse_side=6,forward_sides=[6,6,8],width=.6,material_bounds=[0,8],initial_chi=.25,
        material_values='firstgeometrychi8; secondgeometrychi2and8; analytic volume means at interfaces',
        observations=dict(sources=60,electric_receiver_channels=180,relative_complex_RMS_noise=.01,
                           acquisition='unchanged R8 VectorPR3D three circles of electric dipoles'),
        baseline=dict(cycles=120,material_inner=12,mf='1 first60,2 last60',
                      beta_material='6 first60,1 last60',dual_beta_current='6 first60,3 last60',
                      prior='unchanged lam2e-4 mean neighbour squared difference + l2=1e-6 mean chi²',
                      APCS_coordinates='inherited whiten_apcs=True SVD range of projected IFFT; not literal2019/nativegamma',
                      physical_value_cache='existing exact R8 cache enabled for paid-field reuse'),
        projection=dict(steps=3000,audit_every=50,scale=6.,omega=1.,grouped_ports=True,
                        CPU_cap_seconds=320.,cap_includes_setup=True,
                        ball='one raw iid circular-complex99% full global event; no source whitening'),
        corrector=dict(maxiter=1000,CPU_cap_seconds=60.,ftol=1e-13,gtol=1e-9,maxls=30,
                        scope='same original source-normalized physical Phi; proposal allowances separate'),
        selection='minimum original source-normalized R8 physical Phi including unchanged prior, noGT',
        CPU_prefix_seconds=[10,30,60,120],
        qualification='reduced highindex finiteCM+RR dyadicDDA; coarsecells not continuumaccuracy qualification; notliteral2019reproduction',
        cell_phase_max=dict(inverse6=float(2*np.pi*(.6/6)*np.sqrt(9)),forward8=float(2*np.pi*(.6/8)*np.sqrt(9))),
        seal='all12 finalpredictions plusprefixchoices before offlineGTreads',truth_access='GENERATOR_ONLY_BEFORE_SOLVE')
    dump(HERE/'FROZEN_PROTOCOL.json',protocol)
    (HERE/'data_public').mkdir();(HERE/'sealed_truth').mkdir()
    cpu,wall=time.process_time(),time.perf_counter();parts=geometries(protocol['seed'])
    dump(HERE/'sealed_truth/GEOMETRIES.json',parts)
    inverse=HighContrastVectorPR3D(ns=6);rows=[];clean_by_geometry={}
    for i,cid in enumerate(CASES):
        name='random_high8' if cid.endswith('crossgrid8') else cid
        forward=inverse if not cid.endswith('crossgrid8') else HighContrastVectorPR3D(ns=8)
        chi,integral=volume_average(forward,parts[name]);before=forward.counts.copy()
        clean,J,E=forward.physical(chi)
        forward_counts={k:v-before.get(k,0) for k,v in forward.counts.items()}
        inv_chi,inv_integral=volume_average(inverse,parts[name])
        if cid.endswith('crossgrid8'):
            assert np.array_equal(forward.receivers,inverse.receivers)
            assert np.array_equal(forward.polarizations,inverse.polarizations)
            coarse_clean=clean_by_geometry[name]
        else:
            clean_by_geometry[name]=clean.copy();coarse_clean=clean
        rng=np.random.default_rng(protocol['seed']+100+i)
        sigma=.01*np.linalg.norm(clean)/np.sqrt(clean.size)
        noise=sigma*(rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape))/np.sqrt(2)
        np.savez_compressed(HERE/'data_public'/f'{cid}.npz',Y=clean+noise,sigma=sigma,
                             grid_shape=np.array([6,6,6]),material_bounds=np.array([0.,8.]))
        np.savez_compressed(HERE/'sealed_truth'/f'{cid}.npz',chi=inv_chi,generation_chi=chi,
                             generation_side=forward.ns,Y_clean=clean,Y_coarse_clean=coarse_clean,noise=noise)
        rows.append(dict(case_id=cid,generation_side=forward.ns,inverse_side=6,
            independent_geometry=name,sigma=float(sigma),observation_shape=list(clean.shape),
            forward_counts=forward_counts,forward_setup=forward.setup,
            generation_integral=integral,inverse_integral=inv_integral,
            truth_access='GENERATOR_ONLY'))
    dump(HERE/'GENERATION_RECEIPT.json',dict(status='COMPLETE',rows=rows,
        shared_inverse_model_setup=inverse.setup,cpu_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-wall,qualification=protocol['qualification']))
    print(json.dumps(dict(status='GENERATED_TRUE3D',cases=CASES,cpu=time.process_time()-cpu)),flush=True)


if __name__=='__main__':main()
