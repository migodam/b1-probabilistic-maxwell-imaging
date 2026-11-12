"""Copy actual Austria public observations and enrich acquisition metadata.

No material GT load, new target generation, physical solve or new validation.
Existing truth archives are copied opaquely to separate offline-only storage.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time,shutil,numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.physics_inputs import existing_r8_root
from research.research_max_struct_r8.src.paper_physics import acquisition
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump

def main():
    start,cpu=time.perf_counter(),time.process_time();out=HERE/'classic_public';out.mkdir(exist_ok=False)
    truthout=HERE/'classic_sealed_truth';truthout.mkdir(exist_ok=False)
    r8=existing_r8_root();r12=REPO/'research/research_max_struct_r12'
    contracts=json.loads((r8/'delegated/paper_baselines/PAPER_CONTRACTS.json').read_text())
    caseid='austria_eps4_400MHz';index=next(i for i,c in enumerate(contracts['NIE2016']['cases']) if c['id']==caseid)
    sources=[dict(id='austria_150to32',path=r8/f'RESULTS/paper_benchmarks/{caseid}/PUBLIC_DATA.npz',
        receipt=r8/f'RESULTS/paper_benchmarks/{caseid}/DATA_RECEIPT.json',
        truth=r8/f'data/paper_benchmarks_private/{caseid}.npz',forward=150,maxiter=100),
        dict(id='austria_matched32',path=r12/'delegated/external_matched32/data_public/OBSERVATIONS.npz',
        receipt=r12/'delegated/external_matched32/GENERATION_RECEIPT.json',
        truth=r12/'delegated/external_matched32/sealed_truth/GT.npz',forward=32,maxiter=300)]
    public=acquisition(32,400e6);rows=[]
    for row in sources:
        with np.load(row['path'],allow_pickle=False) as z:Y=z['Y'].copy()
        if Y.shape!=(40,20):raise ValueError('actual public observation dimensions changed')
        old=json.loads(row['receipt'].read_text())
        if row['forward']==150:
            # Known source generator direction and normalization permit the
            # scalar noise norm to be recovered from public Y alone. No clean
            # observation or GT is created/read/saved.
            rng=np.random.default_rng(8600+index)
            direction=rng.normal(size=Y.shape)+1j*rng.normal(size=Y.shape)
            direction/=np.linalg.norm(direction)
            norm2=float(np.vdot(Y,Y).real);c=float(np.vdot(direction,Y).real)
            noise_norm=(-.02*c+np.sqrt((.02*c)**2+4*.99*.01*norm2))/(2*.99)
            sigma_source='known original seed8600+caseindex and fixed-norm source equation; recovered scalar norm from publicY only'
        else:
            noise_norm=float(old['coarse_GT_observed_relative_residual'])*np.linalg.norm(Y)
            sigma_source='existing matched-generation receipt relative_noise_to_noisyY; no cleanY or GT read'
        sigma=float(noise_norm/np.sqrt(Y.size));delta=float(.12*np.linalg.norm(Y))
        arrays=dict(G=public['G'],E0=public['E0'],S=public['S'],grid_shape=np.array([32,32]),
            coords=public['centers'],centers=public['centers'],radius=np.array(public['radius']),
            cell_radius=np.array(public['radius']),k=np.array(public['k']),components=np.array(1),
            kind=np.array('classic'),Y=Y,sigma=np.array(sigma),delta=np.array(delta),
            lower=np.array(0.),upper=np.array(20.),m_init=np.zeros(1024),
            frequency=np.array(400e6),forward_side=np.array(row['forward']),inverse_side=np.array(32),
            baseline_dpcs_rank=np.array(15),baseline_mf=np.array([5,8,10]),
            baseline_beta=np.array([8.,2.,.5]),baseline_beta_current=np.array([8.,8.,8.]),
            baseline_maxiter=np.full(3,row['maxiter']),baseline_coefficient_tol=np.array([1e-4,1e-6,1e-6]),
            noise_sigma_known=np.array(False),
            sigma_role=np.array('RMS of fixed-norm complex noise direction; not iid CN standard deviation'),
            delta_role=np.array('inherited .12 Frobenius discrepancy radius; not iidCN99 certificate'))
        np.savez_compressed(out/f"{row['id']}.npz",**arrays)
        shutil.copyfile(row['path'],out/f"{row['id']}_ORIGINAL_PUBLIC_OBSERVATIONS.npz")
        shutil.copyfile(row['receipt'],out/f"{row['id']}_SOURCE_RECEIPT.json")
        if not row['truth'].is_file():raise FileNotFoundError(row['truth'])
        shutil.copyfile(row['truth'],truthout/f"{row['id']}_OFFLINE_ONLY.npz")
        metadata=dict(case_id=row['id'],source_public=str(row['path']),source_receipt=str(row['receipt']),
            output=str(out/f"{row['id']}.npz"),forward_side=row['forward'],inverse_side=32,
            physical_contract='2m square,20 uniform plane waves,40 receivers at3.75m,400MHz,epsilon_r4/chi3',
            quadrature='existing equal-area disk pulse h/sqrt(pi),analytic self; not exact square integration',
            initialization='all1024unknown cells exactly real chi0, original uniform air host',
            material_box=[0,20],noise_relative_to_clean=.1,noise_fixed_norm=True,
            actual_noise_RMS_per_channel=sigma,noise_norm_reconstruction=sigma_source,
            delta=delta,delta_relative_to_observed=.12,
            statistical_noise_certificate=False,
            prior=dict(smooth=1e-4,tv=1e-4,epsilon=.01),
            baseline=dict(dpcs_rank=15,mf=[5,8,10],beta_material=[8,2,.5],
                beta_current_DualGlobal=[8,8,8],beta_current_CIE=[8,2,.5],
                maxiter_per_stage=row['maxiter'],coefficient_tol=[1e-4,1e-6,1e-6],
                sourcewise_CIE_Dual=True,global_PR_control=True,initial_APIC='zero native Fourier coefficients',
                material_inner_maxiter=60,physical_warm_start_historical=False),
            original_NIE2016_contract=dict(source=str(r8/'delegated/paper_baselines/PAPER_CONTRACTS.json'),
                DOI=contracts['NIE2016']['doi'],forward_side=150,inverse_side=64,
                additional_prior=False,initial_contrast=0,initial_APIC=0,
                stop_delta2D=[.01,.001,.001]),
            adaptation='inverse32 instead of paper64; disk quadrature; real finitebox; common prior/corrector; not author-code reproduction',
            scientific_scope='reused crossgrid classic challenge' if row['forward']==150 else 'post-exposure matched-grid diagnostic, not fresh confirmation',
            truth_copy=dict(source=str(row['truth']),destination=str(truthout/f"{row['id']}_OFFLINE_ONLY.npz"),
                method='opaque byte copy; no truth array opened; offline scoring only after all prediction/choice seals'),
            new_physics=dict(forward=0,LU=0,RHS=0,validation=0,Green_assembly=1 if row is sources[0] else 0))
        dump(out/f"{row['id']}_METADATA.json",metadata);rows.append(metadata)
    dump(out/'PACKAGING_RECEIPT.json',dict(status='COMPLETE',rows=rows,new_physical_solves=0,
        new_LU=0,new_RHS=0,new_validation=0,acquisition_assembly_count=1,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,
        truth_access='opaque file copy only; no truth values loaded',historical_results_modified=False))
    print('Two actual Austria public inputs copied and enriched; zero new physical solves',flush=True)
if __name__=='__main__':main()
