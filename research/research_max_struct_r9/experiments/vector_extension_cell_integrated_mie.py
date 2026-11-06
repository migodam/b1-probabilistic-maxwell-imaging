"""One fixed IGT Mie16 paired control, public sphere qualification only.
No inverse GT, material/frequency fitting, old-backend edits or GPU use.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import time
WALL=time.perf_counter()
from pathlib import Path
import sys,json,argparse,traceback,datetime,importlib.util
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
CONFIG=ROOT/'configs/vector_extension_cell_integrated_mie16.json';OUT=ROOT/'RESULTS/vector_extension/cell_integrated_mie16'

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def load(p):return json.loads(Path(p).read_text())
def dump(p,r):Path(p).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')

def main():
    import psutil
    p=argparse.ArgumentParser();p.add_argument('--cpu',type=int,default=8);args=p.parse_args();process=psutil.Process()
    if os.name=='nt':
        process.cpu_affinity([args.cpu]);spec=importlib.util.spec_from_file_location('r9_windows_job',ROOT/'delegated/3d_extension/windows_resolved.py')
        joblib=importlib.util.module_from_spec(spec);spec.loader.exec_module(joblib);job=joblib.job_object()
    else:job=dict(active=False,reason='nonWindowsdiagnostic')
    import numpy as np
    from scipy import linalg as la
    from research.research_max_struct_r9.src.cell_integrated_green import CellIntegratedGreen,self_scalar,integrated_blocks,cube_rule
    from research.research_max_struct_r1.src.vector_backend import plane_waves
    cfg=load(CONFIG);OUT.mkdir(parents=True,exist_ok=False);dump(OUT/'FROZEN_CONFIG.json',cfg)
    r=dict(status='RUNNING',started_utc=utc(),new_inverse_runs=0,truth_access='NO_INVERSE_GT; declaredpublicspherebenchmarkonly',GPU=False,
        platform='XINAN_Windows_i9' if os.name=='nt' else sys.platform,CPU_affinity=process.cpu_affinity() if hasattr(process,'cpu_affinity') else None,
        job_object=job,one_BLAS_thread_environment={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')},
        wholeprocess_CPU_source='psutiluser+system/GetProcessTimes onWindows',counts={},primary_source='2019 CIE Eq2 volumeintegral; rootderived cubicself trace/delta',
        interpretation='fixed standardpulsecollocationVIE diagnostic, not a new theorem, passivitycertificate or fitted normalization')
    dump(OUT/'RECEIPT.json',r)
    try:
        with np.load(ROOT/cfg['reference'],allow_pickle=False) as z:
            points=z['points'].copy();fraction=z['fraction'].copy();x=z['chi'].copy();theta=z['theta'].copy();oldfar=z['far'].copy();reference=z['Mie_far'].copy()
        old=load(ROOT/cfg['old_receipt']);h=cfg['h'];k=cfg['k'];a=h/2
        model=CellIntegratedGreen(h,k,cfg['near_order'],cfg['far_order'],cfg['self_order'],cfg['near_Chebyshev_radius'])
        # One fixed independent regular-part check, not a quadrature-order scan.
        qq,ww=cube_rule(h,cfg['checks']['self_independent_cartesian_regularpart_order']);rr=np.linalg.norm(qq,axis=1)
        static_integral=a*a*(24*np.log((1+np.sqrt(3))/np.sqrt(2))-2*np.pi)
        independent_integral=static_integral+np.sum(ww*np.expm1(1j*k*rr)/rr)
        independent_self=-1/3+k*k/(6*np.pi)*independent_integral
        self20=self_scalar(h,k,cfg['checks']['self_higher_Duffy_order'])
        near=np.asarray(cfg['checks']['near_fixed_offset']);far=np.asarray(cfg['checks']['far_fixed_offset'])
        near8=integrated_blocks(near[None]*h,h,k,cfg['near_order'])[0];near12=integrated_blocks(near[None]*h,h,k,cfg['checks']['near_reference_order'])[0]
        far3=integrated_blocks(far[None]*h,h,k,cfg['far_order'])[0];far5=integrated_blocks(far[None]*h,h,k,cfg['checks']['far_reference_order'])[0]
        v=np.asarray(cfg['checks']['symmetry_fixed_offset']);perm=np.array([1,0,2]);sym=integrated_blocks(np.array([v,-v,v[perm]])*h,h,k,cfg['near_order'])
        checks=dict(self_primary=[model.self_density.real,model.self_density.imag],self_Duffy20=[self20.real,self20.imag],self_independent_cartesian_staticplusregular=[independent_self.real,independent_self.imag],
            self16_vs20_absolute=abs(model.self_density-self20),self16_vs_independent_absolute=abs(model.self_density-independent_self),
            static_cube_Coulomb_integral=static_integral,static_integral_formula='a²[24log((1+sqrt3)/sqrt2)-2pi]',
            near8_vs12_relative=float(la.norm(near8-near12)/la.norm(near12)),far3_vs5_relative=float(la.norm(far3-far5)/la.norm(far5)),
            inversion_symmetry_relative=float(la.norm(sym[0]-sym[1])/la.norm(sym[0])),axis_permutation_relative=float(la.norm(sym[2]-sym[0][perm][:,perm])/la.norm(sym[0])),
            block_reciprocity_relative=float(la.norm(sym[0]-sym[0].T)/la.norm(sym[0])),scope='threefixedspotorders plusonesymmetrylag; no resultdrivenorderchange')
        dump(OUT/'QUADRATURE_CHECK.json',checks)
        coordinates=np.rint((points-points[0])/h).astype(int)
        if np.max(abs((points-points[0])/h-coordinates))>1e-12:raise ValueError('Saved benchmark is not same uniform grid')
        G=model.assemble_dense(coordinates);d=np.repeat(x,3);inc=plane_waves(points,k,[[0,0,1]],[[1,0,0]])
        K=np.eye(len(d),dtype=complex)-d[:,None]*G;factor=la.lu_factor(K,check_finite=False);r['counts']['physical_LU']=1
        J=la.lu_solve(factor,d[:,None]*inc,check_finite=False);r['counts']['physical_source_RHS']=1
        E=inc+G@J;r['counts']['physical_G_columns']=1
        directions=np.c_[np.sin(theta),theta*0,np.cos(theta)];newfar=model.far_field(points,J,directions)[:,:,0]
        error=float(la.norm(newfar-reference)/la.norm(reference));old_error=float(la.norm(oldfar-reference)/la.norm(reference))
        mu,w=np.polynomial.legendre.leggauss(24);phi=np.arange(32)*2*np.pi/32;mm,pp=np.meshgrid(mu,phi,indexing='ij')
        directions_quad=np.stack([np.sqrt(1-mm*mm)*np.cos(pp),np.sqrt(1-mm*mm)*np.sin(pp),mm],axis=-1).reshape(-1,3)
        q=model.far_field(points,J,directions_quad)[:,:,0].reshape(24,32,3);cross=float(np.sum(w[:,None]*np.sum(abs(q)**2,axis=2))*2*np.pi/32)
        sigma=old['mie_scattering_cross_section'];gain=np.vdot(reference,newfar)/np.vdot(reference,reference)
        r.update(status='COMPLETED',nside=16,active_cells=len(points),vector_dimension=len(d),epsilon_r=9,chi=8,radius=.145,h=h,
            geometry_volume=float(fraction.sum()*h**3),geometry_fraction_rule='exactsame savedMG24fractions andchi asold16',
            integrated_far_complex_vector_relative_error=error,retained_old_CMRR_far_error=old_error,old_CMRR_rerun=False,
            scattering_cross_section=cross,mie_scattering_cross_section=sigma,scattering_cross_section_relative_error=abs(cross-sigma)/sigma,
            retained_old_CMRR_cross_section_error=old['scattering_cross_section_relative_error'],
            linear_relative_residual=float(la.norm(J-d[:,None]*E)/la.norm(d[:,None]*inc)),
            far_cube_form_factor_applied=True,far_cube_form_factor_formula='product np.sinc(k*h*n/(2*pi)); h³current integral',
            unapplied_diagnostic_gain=[float(gain.real),float(gain.imag)],unapplied_gain_shape_residual=float(la.norm(newfar-gain*reference)/la.norm(reference)),
            operator=model.receipt,quadrature_checks=checks,continuum_certificate='NOT_CLAIMED; one pairedresolutionpoint, sameapproximateMGboundary')
        np.savez_compressed(OUT/'fields.npz',points=points,fraction=fraction,chi=x,J=J,E=E,theta=theta,far=newfar,Mie_far=reference,retained_old_CMRR_far=oldfar)
    except Exception as exc:r.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
    cp=process.cpu_times();mem=process.memory_info();r.update(finished_utc=utc(),wholeprocess_CPU_seconds=cp.user+cp.system,script_wall_seconds=time.perf_counter()-WALL,
        process_peak_RSS_bytes=getattr(mem,'peak_wset',None),peak_memory_source='Windowspeak_wset' if hasattr(mem,'peak_wset') else 'UNKNOWN_NOT_ZERO')
    dump(OUT/'RECEIPT.json',r);print(json.dumps(r),flush=True);return 0 if r['status']=='COMPLETED' else 1

if __name__=='__main__':raise SystemExit(main())
