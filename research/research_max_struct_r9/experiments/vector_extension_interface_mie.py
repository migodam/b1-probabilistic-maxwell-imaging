"""Two fixed known-interface Mie diagnostics, two LU/RHS; no inverse GT.
Workers do not load Mie outputs. Parent seals both before one offline score.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import time
WALL=time.perf_counter()
from pathlib import Path
import sys,json,argparse,datetime,subprocess,traceback,importlib.util
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
CONFIG=ROOT/'configs/vector_extension_interface_mie16.json';OUT=ROOT/'RESULTS/vector_extension/interface_mie16'
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def load(p):return json.loads(Path(p).read_text())
def dump(p,r):Path(p).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def stamp(p):
    p=Path(p);s=p.stat();return dict(path=str(p.relative_to(REPO)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def cpu(process):
    c=process.cpu_times();return c.user+c.system

def worker(cfg,arm,process):
    import numpy as np
    from scipy import linalg as la
    from research.research_max_struct_r9.src.cell_integrated_green import CellIntegratedGreen
    from research.research_max_struct_r1.src.vector_backend import plane_waves
    dest=OUT/arm;dest.mkdir(parents=True,exist_ok=False)
    r=dict(status='RUNNING',arm=arm,classification=cfg['classification'],started_utc=utc(),counts={},inverse_GT_access=False,
        Mie_output_access='NONE; no score/reference values loaded byworker',platform='XINAN_Windows_i9',CPU_affinity=process.cpu_affinity(),GPU=False)
    dump(dest/'RECEIPT.json',r)
    try:
        # NPZ lazy key access deliberately excludes Mie_far, far and old chi.
        with np.load(ROOT/cfg['geometry_source'],allow_pickle=False) as z:points=z['points'].copy();f=z['fraction'].copy()
        N=len(points);M=3*N;h=cfg['h'];k=cfg['k'];eps0=cfg['epsilon_material'];eye=np.eye(3)
        et=1+f*(eps0-1);epsilon=et[:,None,None]*eye
        if arm=='normal_tensor':
            cut=(f>0)&(f<1);norm=np.linalg.norm(points[cut],axis=1)
            if np.any(norm==0):raise ValueError('Undefined cutcell sphere normal')
            n=points[cut]/norm[:,None];en=1/(1-f[cut]+f[cut]/eps0)
            epsilon[cut]+=(en-et[cut])[:,None,None]*n[:,:,None]*n[:,None,:]
            epsilon[f==1]=eps0*eye
        elif arm!='scalar_arithmetic':raise ValueError('Unknown frozen arm')
        chi=epsilon-eye;kernel=CellIntegratedGreen(h,k,**dict(near_order=8,far_order=3,self_order=16,near_radius=2))
        coordinates=np.rint((points-points[0])/h).astype(int);G=kernel.assemble_dense(coordinates)
        inc=plane_waves(points,k,[[0,0,1]],[[1,0,0]])
        # Exact R1 TensorDDA E-form L=I-G*blockdiag(chi), same dense LU;
        # block column multiplication avoids another dense blockdiag GEMM.
        L=-np.einsum('iajb,jbc->iajc',G.reshape(N,3,N,3),chi,optimize=False).reshape(M,M)
        L.flat[::M+1]+=1;factor=la.lu_factor(L,check_finite=False);r['counts']['physical_LU']=1
        E=la.lu_solve(factor,inc,check_finite=False);r['counts']['physical_source_RHS']=1
        J=np.einsum('nab,nbp->nap',chi,E.reshape(N,3,-1)).reshape(M,-1)
        field=inc+G@J;r['counts']['physical_G_columns']=1
        theta=np.linspace(0,np.pi,181);directions=np.c_[np.sin(theta),theta*0,np.cos(theta)]
        far=kernel.far_field(points,J,directions)[:,:,0]
        mu,w=np.polynomial.legendre.leggauss(24);phi=np.arange(32)*2*np.pi/32;mm,pp=np.meshgrid(mu,phi,indexing='ij')
        quad=np.stack([np.sqrt(1-mm*mm)*np.cos(pp),np.sqrt(1-mm*mm)*np.sin(pp),mm],axis=-1).reshape(-1,3)
        far_quad=kernel.far_field(points,J,quad)[:,:,0].reshape(24,32,3)
        np.savez_compressed(dest/'PREDICTION.npz',points=points,fraction=f,epsilon=epsilon,E=E,J=J,theta=theta,far=far,far_quad=far_quad,quad_weights=w)
        r.update(status='SEALED_UNSCORED_OFFLINE_DIAGNOSTIC',active_cells=N,vector_dimension=M,operator=kernel.receipt,
            formulation='R1 TensorDDA tensor E-form, I-G*chi; block3x3 assembly and inheritedSciPydenseLU',
            linear_relative_residual=float(la.norm(E-field)/la.norm(inc)),far_cube_sinc_applied=True)
    except Exception as exc:r.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
    r.update(finished_utc=utc(),wholeprocess_CPU_seconds=cpu(process),script_wall_seconds=time.perf_counter()-WALL,peak_RSS_bytes=process.memory_info().peak_wset)
    dump(dest/'RECEIPT.json',r);print(json.dumps(r),flush=True);return 0 if r['status']=='SEALED_UNSCORED_OFFLINE_DIAGNOSTIC' else 1

def seal_score(cfg):
    import numpy as np
    from scipy import linalg as la
    entries=[]
    for arm in cfg['arms']:
        folder=OUT/arm;r=load(folder/'RECEIPT.json')
        if r['status']!='SEALED_UNSCORED_OFFLINE_DIAGNOSTIC':raise RuntimeError('Both terminal predictions required before Mie score')
        entries.append(dict(arm=arm,files=[stamp(folder/'PREDICTION.npz'),stamp(folder/'RECEIPT.json')]))
    seal=dict(status='BOTH_PREDICTIONS_SEALED_BEFORE_OFFLINE_MIE_SCORING',created_utc=utc(),entries=entries,evidence=[stamp(CONFIG),stamp(Path(__file__)),stamp(ROOT/'src/cell_integrated_green.py')],Mie_reference_access='NONE_BEFORE_THIS_SEAL',inverse_GT_access=False)
    dump(OUT/'PREDICTIONS_SEAL.json',seal)
    # First Mie output access in this newtwo-arm experiment, after both seals.
    with np.load(ROOT/cfg['geometry_source'],allow_pickle=False) as z:reference=z['Mie_far'].copy()
    old=load(ROOT/'RESULTS/vector_extension/validation_refinement16/RECEIPT.json');scores=[]
    for arm in cfg['arms']:
        with np.load(OUT/arm/'PREDICTION.npz',allow_pickle=False) as z:
            far=z['far'];q=z['far_quad'];w=z['quad_weights'];eps=z['epsilon'];f=z['fraction']
        eig=np.linalg.eigvalsh(eps);cross=float(np.sum(w[:,None]*np.sum(abs(q)**2,axis=2))*2*np.pi/32)
        scores.append(dict(arm=arm,far_complex_relative_L2_error=float(la.norm(far-reference)/la.norm(reference)),scattering_cross_section=cross,
            scattering_cross_section_relative_error=abs(cross-old['mie_scattering_cross_section'])/old['mie_scattering_cross_section'],
            eigenvalue_min=float(eig.min()),eigenvalue_max=float(eig.max()),real_SPD=bool(np.min(eig)>0),
            reciprocal_symmetry_Linf=float(np.max(abs(eps-eps.transpose(0,2,1)))),fullcell_epsilon9_max_error=float(np.max(abs(eps[f==1]-9*np.eye(3)))),
            cutcells=int(np.sum((f>0)&(f<1))),Mie_assessment='ONE_OFFLINE_PASS_AFTER_BOTH_SEALED',inverse_recovery_claim=False))
    dump(OUT/'SCORES.json',dict(classification=cfg['classification'],scores=scores,old_MG_IGT_error=0.4339336694219462,
        no_old_prediction_or_score_changes=True,scope='knownsphereinterface materialrepresentation only; no unknowninterface recovery or uniquecause proof'))

def main():
    import psutil
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker']);p.add_argument('--arm');args=p.parse_args();process=psutil.Process();process.cpu_affinity([8]);cfg=load(CONFIG)
    if args.mode=='worker':return worker(cfg,args.arm,process)
    spec=importlib.util.spec_from_file_location('joblib',ROOT/'delegated/3d_extension/windows_resolved.py');joblib=importlib.util.module_from_spec(spec);spec.loader.exec_module(joblib);job=joblib.job_object()
    OUT.mkdir(parents=True,exist_ok=False);dump(OUT/'FROZEN_CONFIG.json',cfg);jobs=[]
    for arm in cfg['arms']:
        command=[sys.executable,str(Path(__file__).resolve()),'worker','--arm',arm];t=time.perf_counter();started=utc()
        with (OUT/(arm+'.log')).open('x') as log:child=subprocess.Popen(command,cwd=REPO,stdout=log,stderr=subprocess.STDOUT);code=child.wait()
        jobs.append(dict(arm=arm,command=command,started_utc=started,finished_utc=utc(),exit_code=code,wall_seconds=time.perf_counter()-t))
        if code:dump(OUT/'PARENT_RECEIPT.json',dict(status='FAILED',jobs=jobs,parent_CPU_seconds=cpu(process),job_object=job));return code
    seal_score(cfg)
    dump(OUT/'PARENT_RECEIPT.json',dict(status='COMPLETED_OFFLINE_DIAGNOSTIC',jobs=jobs,parent_CPU_seconds=cpu(process),parent_wall_seconds=time.perf_counter()-WALL,
        CPU_affinity=[8],GPU=False,job_object=job,CPU_scope='Parentonlyimports/dispatch/seal/offlineMieassessment; workers separatelycharged',new_inverse_runs=0))
    print(json.dumps(load(OUT/'SCORES.json')),flush=True);return 0

if __name__=='__main__':raise SystemExit(main())
