"""R9 Phase2 vector extension: frozen preparation, public-only workers, seal/score.
Never edits R8/PhaseI. CLI: all | prepare | worker | score | report.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import time
ENTRY_WALL,ENTRY_CPU=time.perf_counter(),time.process_time()
from pathlib import Path
import sys,json,argparse,traceback,subprocess,resource,gc,shutil
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
from research.research_max_struct_r9.src.vector_extension import StrictVectorInverse,ARMS,target_material,paid_density_state
from research.research_max_struct_r1.src.vector_backend import TensorDDA,plane_waves
from research.research_max_struct_r1.experiments.qualify_backend import mie
from collections import Counter
CONFIG=ROOT/'configs/vector_extension_phase2.json';PUBLIC=ROOT/'data_public/vector_extension';TRUTH=ROOT/'sealed_truth/vector_extension';OUT=ROOT/'RESULTS/vector_extension';LOGS=ROOT/'delegated/3d_extension'

def utc():return datetime.now(timezone.utc).isoformat()
def load(p):return json.loads(Path(p).read_text())
def dump(p,x):
    def conv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,Path):return str(v)
        raise TypeError(type(v).__name__)
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(x,indent=2,allow_nan=False,default=conv)+'\n');tmp.replace(p)
def stamp(p):
    p=Path(p);s=p.stat();return dict(path=str(p.relative_to(REPO)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def save_state(p,state):
    values={k:state[k] for k in ('m','J','E','gradient') if k in state}
    values['Y_prediction']=state['Y'];np.savez_compressed(p,**values)
def fieldmeta(model):
    return dict(nside=model.nside,material_dof=model.N,vector_dimension=model.M,domain_width=model.width,h=model.h,
        k=2*np.pi,cell_volume=model.h**3,sources=model.receivers.tolist(),receivers=model.receivers.tolist(),
        source_polarizations=model.polarizations.tolist(),source_circle_radius=model.width,
        source_centers_scale=model.width/3,receiver_components=3,complex_observations=10800,
        epsilon_background=1,mu_r=1,time_convention='exp(-i omega t)',self_term_real=float(model.s.real),self_term_imag=float(model.s.imag))

def qualify(cfg):
    folder=OUT/'validation';folder.mkdir(exist_ok=True)
    if (folder/'RECEIPT.json').exists():raise FileExistsError('Validation receipt already exists')
    wall,cpu=time.perf_counter(),time.process_time();rows=[];failures=[]
    theta=np.linspace(0,np.pi,181);directions=np.c_[np.sin(theta),theta*0,np.cos(theta)];k=2*np.pi
    epsilon=cfg['validation']['epsilon_r'];contrast=epsilon-1;radius=cfg['validation']['radius']
    _,s2,_=mie(radius,k,epsilon,np.cos(theta));reference=1j/k*s2[:,None]*np.c_[np.cos(theta),theta*0,-np.sin(theta)]
    a=(np.arange(24)+.5)/24-.5;sub=np.stack(np.meshgrid(a,a,a,indexing='ij'),-1).reshape(-1,3)
    for nside in cfg['validation']['nside']:
        start,cp=time.perf_counter(),time.process_time();cost=Counter()
        try:
            width=cfg['validation']['domain_width'];h=width/nside;ijk=np.indices((nside,)*3).reshape(3,-1).T;points=(ijk-(nside-1)/2)*h
            backend=TensorDDA(points,h**3);fraction=np.array([np.mean(np.sum((p+h*sub)**2,axis=1)<=radius**2) for p in points])
            x=3*fraction*contrast/(3+(1-fraction)*contrast)
            inc=plane_waves(points,k,[[0,0,1]],[[1,0,0]])
            state=paid_density_state(backend.G,inc,x,None,cost)
            phase=np.exp(-1j*k*directions@points.T);dip=phase@(h**3*state['J'].reshape(-1,3))
            far=k*k/(4*np.pi)*(dip-directions*np.sum(directions*dip,axis=1)[:,None])
            error=float(la.norm(far-reference)/la.norm(reference))
            row=dict(status='COMPLETED',nside=nside,N=nside**3,vector_dimension=3*nside**3,h=h,
                epsilon_r=epsilon,chi=contrast,radius=radius,relative_Mie_far_vector_error=error,
                geometry_volume=float(fraction.sum()*h**3),exact_volume=4*np.pi*radius**3/3,
                linear_relative_residual=state['linear_relative_residual'],counts=dict(cost),
                CPU_seconds=time.process_time()-cp,wall_seconds=time.perf_counter()-start)
            np.savez_compressed(folder/f'grid{nside}.npz',points=points,fraction=fraction,chi=x,far=far,Mie_far=reference,theta=theta)
            rows.append(row);del backend,state,points,far;gc.collect()
        except Exception as exc:
            failure=dict(nside=nside,error=repr(exc),traceback=traceback.format_exc(),counts=dict(cost),CPU_seconds=time.process_time()-cp,wall_seconds=time.perf_counter()-start);failures.append(failure)
    result=dict(status='COMPLETED_WITH_FAILURES' if failures else 'COMPLETED',rows=rows,failures=failures,
        rule='fixed24³midpoints and MGstaticpolarizabilitymix; nofit/rescale; planeincidenceE=x,k=z',
        error_decreased=(rows[-1]['relative_Mie_far_vector_error']<rows[0]['relative_Mie_far_vector_error']) if len(rows)==2 else None,
        continuum_accuracy='NOT_CERTIFIED; measured bounded geometry/contrast qualification, not exactcutcell',
        no_halt_other_inversions=True,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        process_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dump(folder/'RECEIPT.json',result);print(json.dumps(dict(event='MIE_QUALIFICATION',status=result['status'],rows=rows)),flush=True)

def prepare(cfg):
    if (PUBLIC/'DATASET_COMPLETE.json').exists():raise FileExistsError('Dataset already frozen')
    qualify(cfg);wall,cpu=time.perf_counter(),time.process_time();contract=load(TRUTH/'GENERATION_CONTRACT.json')['cases'];rows=[]
    legacy=REPO.parent/'r8_worktree/research/research_max_struct_r8/RESULTS/vector3d/refinement8'
    for spec in cfg['cases']:
        case=spec['case_id'];dest=PUBLIC/case;private=TRUTH/case;dest.mkdir(exist_ok=False);private.mkdir(exist_ok=False)
        t,c=time.perf_counter(),time.process_time()
        if case=='four_cubes_legacy':
            with np.load(legacy/'public_observations.npz',allow_pickle=False) as z:Y=z['Y'].copy();clean=z['clean_Y'].copy()
            with np.load(legacy/'truth_for_scoring.npz',allow_pickle=False) as z:truth=z['chi'].copy()
            model=StrictVectorInverse(8,.6,cfg,auxiliary=False);geometry=fieldmeta(model)
            meta=dict(case_id=case,noise_model='legacy_fixed_energy_bounded',noise_radius=float(.1/.9*la.norm(Y)),
                radius_scope='bounds original fixedenergy noise only; NOT fine/coarsemodeldifference',acquisition=geometry,
                generation='copied R8publicY unchanged; no repeatedfine10forward',phase2_new_simulation_counts=dict(physical_LU=0,physical_source_RHS=0))
            simulation=load(legacy/'SIMULATION.json');mismatch=simulation['crossgrid_forward_discrepancy']
            counts=dict(physical_LU=0,physical_source_RHS=0);setup=model.setup;del model;gc.collect()
        else:
            fine=StrictVectorInverse(10,spec['domain_width'],cfg,auxiliary=False)
            xfine=target_material(case,fine.points,fine.h,contract);st=fine.state(xfine);clean=st['Y'].copy()
            fine_counts=dict(fine.counts);fine_setup=fine.setup;del fine,st,xfine;gc.collect()
            coarse=StrictVectorInverse(8,spec['domain_width'],cfg,auxiliary=False)
            truth=target_material(case,coarse.points,coarse.h,contract);ref=coarse.state(truth)
            mismatch=float(la.norm(ref['Y']-clean)/la.norm(clean));geometry=fieldmeta(coarse)
            sigma=float(spec['expected_relative_noise']*la.norm(clean)/np.sqrt(clean.size))
            rng=np.random.default_rng(contract[case]['noise_seed']);noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape));Y=clean+noise
            meta=dict(case_id=case,noise_model='iid_circular_complex',sigma=sigma,
                sigma_convention='E|n|²=sigma² commonrawcomplexobservation',expected_relative_Frobenius_noise=.03,
                acquisition=geometry,
                noise_rule_scope='CNnoiseball only, NOT discretizationerrorbudget')
            realization_noise_ratio=float(la.norm(noise)/la.norm(clean))
            counts=dict(fine=fine_counts,coarse=dict(coarse.counts));setup=dict(fine=fine_setup,coarse=coarse.setup)
            del coarse,ref,noise;gc.collect()
        # Public payload deliberately contains no clean/noise/material/support/seed.
        np.savez_compressed(dest/'observations.npz',Y=Y);dump(dest/'PUBLIC.json',meta)
        np.savez_compressed(private/'truth.npz',chi=truth,Y_clean=clean)
        row=dict(case_id=case,status='GENERATED_OR_COPIED',crossgrid_clean_relative_difference=mismatch,
            mismatch_scope='independentfine10/coarse8assembly at declaredcoarsereference; sameDDAfamily, notindependentPDE',
            simulation_counts=counts,operator_setup=setup,actual_new_iid_noise_ratio=(realization_noise_ratio if case!='four_cubes_legacy' else None),CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
        dump(OUT/'generation'/f'{case}.json',row);rows.append(row);del truth,clean,Y
        print(json.dumps(dict(event='DATASET_PREPARED',case=case,mismatch=mismatch)),flush=True)
    dump(PUBLIC/'DATASET_COMPLETE.json',dict(status='PUBLIC_DATASET_FROZEN',created_utc=utc(),cases=[x['case_id'] for x in cfg['cases']],truth_arrays='sealed_truth/vector_extension only',self_contained_public=True))
    dump(OUT/'generation/RECEIPT.json',dict(status='COMPLETE',rows=rows,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,legacy_new_forward_LU=0,truth_access='GENERATOR_ONLY',solver_access_to_truth=False))

def worker(cfg,case,arm):
    out=OUT/case/arm;out.mkdir(parents=True,exist_ok=False)
    receipt=dict(status='RUNNING',case_id=case,arm=arm,started_utc=utc(),truth_access='NONE',scoring='NOT_RUN',
        config=stamp(CONFIG),input=[stamp(PUBLIC/case/'PUBLIC.json'),stamp(PUBLIC/case/'observations.npz')],
        implementation='R9strictcurrent/material orchestration over inherited R1/R8 fullvectorDDA; no newbackend',
        prior_phase_files_modified=False,threads_per_process=1)
    dump(out/'RECEIPT.json',receipt);dump(out/'FROZEN_CONFIG.json',cfg)
    try:
        public=load(PUBLIC/case/'PUBLIC.json')
        with np.load(PUBLIC/case/'observations.npz',allow_pickle=False) as z:
            if z.files!=['Y']:raise ValueError('Unexpected private/publicpayload')
            Y=z['Y'].copy()
        model=StrictVectorInverse(8,public['acquisition']['domain_width'],cfg,auxiliary=arm!='direct')
        def progress(row):print(json.dumps(dict(case_id=case,arm=arm,**row)),flush=True)
        arrays,run=model.run(Y,public,arm,progress)
        save_state(out/'initial.npz',arrays['initial']);save_state(out/'raw.npz',arrays['raw']);save_state(out/'prediction.npz',arrays['final'])
        for name,state in arrays['snapshots'].items():save_state(out/(name+'.npz'),state)
        for st in arrays['stages']:
            np.savez_compressed(out/f"stage_{st['stage']:02d}_joint.npz",m=st['m'],J=st['J'],Q=st['Q'],gamma=st['gamma'])
            dump(out/f"stage_{st['stage']:02d}_joint.json",{k:v for k,v in st.items() if k not in ('m','J','Q','gamma')})
        dump(out/'PAID_BUDGET.json',arrays['budget']);np.savez_compressed(out/'PAID_BUDGET_MATERIAL.npz',m=arrays['budget_material'])
        dump(out/'JOINT_HISTORY.json',run.pop('joint_history'))
        receipt.update(run,operator_setup=model.setup,phase2_incremental_CPU_seconds=time.process_time()-ENTRY_CPU,
            phase2_incremental_wall_seconds=time.perf_counter()-ENTRY_WALL,
            oneoff_setup_and_online_CPU_seconds=model.setup['CPU_seconds']+run['online_CPU_seconds'],
            cost_scope='Wholeworker CPU includes imports/setup/online/export; nestedsetup/online overlapwithtotal; no borrowedR8warmcost')
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),counts=dict(model.counts) if 'model' in locals() else {},
            phase2_incremental_CPU_seconds=time.process_time()-ENTRY_CPU,phase2_incremental_wall_seconds=time.perf_counter()-ENTRY_WALL)
    receipt.update(finished_utc=utc(),process_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dump(out/'RECEIPT.json',receipt);print(json.dumps(dict(event='WORKER_END',case_id=case,arm=arm,status=receipt['status'],error=receipt.get('error'))),flush=True)
    return 0 if receipt['status']=='SEALED_PUBLIC_PREDICTION' else 1

def run_queue(cfg):
    folder=LOGS/'phase2_queue';folder.mkdir(exist_ok=False);pending=[]
    # First launch different geometries/methods; thereafter method/case rotation.
    cases=[c['case_id'] for c in cfg['cases']]
    for i,arm in enumerate(cfg['arms']):
        for j in range(len(cases)):pending.append((cases[(i+j)%len(cases)],arm))
    start=time.perf_counter();active={};jobs=[]
    while pending or active:
        while pending and len(active)<2:
            case,arm=pending.pop(0);path=folder/(case+'__'+arm+'.log');fh=path.open('x');command=[sys.executable,str(Path(__file__).resolve()),'worker','--case',case,'--arm',arm]
            t=time.perf_counter();p=subprocess.Popen(command,cwd=REPO,stdout=fh,stderr=subprocess.STDOUT)
            row=dict(case_id=case,arm=arm,command=command,log=str(path.relative_to(REPO)),started_utc=utc(),pid=p.pid)
            active[p]=(fh,t,row);print(json.dumps(dict(event='START',**row)),flush=True)
        for p,(fh,t,row) in list(active.items()):
            code=p.poll()
            if code is None:continue
            fh.close();row.update(returncode=code,ended_utc=utc(),wall_seconds=time.perf_counter()-t)
            recpath=OUT/row['case_id']/row['arm']/'RECEIPT.json'
            rec=load(recpath) if recpath.exists() else {};row.update(status=rec.get('status','MISSING_RECEIPT'),error=rec.get('error'))
            if code:row['log_tail']=(REPO/row['log']).read_text(errors='replace')[-6000:]
            jobs.append(row);del active[p];print(json.dumps(dict(event='END',**row)),flush=True)
        dump(folder/'RECEIPT.json',dict(status='RUNNING' if active or pending else 'COMPLETE',jobs=jobs,pending=len(pending),running=len(active),max_processes=2,wall_seconds=time.perf_counter()-start,truth_access='NONE',scoring='NOT_RUN'))
        if active:time.sleep(.5)
    return jobs

def seal_and_score(cfg):
    sealpath=OUT/'PUBLIC_PREDICTIONS_SEAL.json'
    if (OUT/'SCORES.json').exists():raise FileExistsError('Phase2 scores already exist')
    entries=[]
    for case in [s['case_id'] for s in cfg['cases']]:
        for arm in cfg['arms']:
            d=OUT/case/arm;r=load(d/'RECEIPT.json')
            if r['status'] not in ('SEALED_PUBLIC_PREDICTION','FAILED'):raise RuntimeError('Nonterminal prediction')
            if r['status']=='SEALED_PUBLIC_PREDICTION' and not (d/'prediction.npz').exists():raise RuntimeError('Missingprediction')
            entries.append(dict(case_id=case,arm=arm,status=r['status'],files=[stamp(f) for f in sorted(d.iterdir()) if f.is_file()]))
    evidence=[stamp(CONFIG),stamp(PUBLIC/'DATASET_COMPLETE.json')]
    for case in [s['case_id'] for s in cfg['cases']]:evidence.extend([stamp(PUBLIC/case/'PUBLIC.json'),stamp(PUBLIC/case/'observations.npz')])
    seal=dict(status='ALL24_PUBLIC_PREDICTIONS_SEALED',created_utc=utc(),entries=entries,evidence=evidence,truth_access='NONE',hash_checks='NONE')
    if sealpath.exists():raise FileExistsError(sealpath)
    dump(sealpath,seal)
    for s in evidence+[v for e in entries for v in e['files']]:
        if stamp(REPO/s['path'])!=s:raise RuntimeError('Changedpublicevidence')
    # First private scoring reads occur only after all24 terminal public outputs.
    scores=[]
    for case in [s['case_id'] for s in cfg['cases']]:
        with np.load(TRUTH/case/'truth.npz',allow_pickle=False) as z:truth=z['chi'].copy()
        norm=la.norm(truth)
        for arm in cfg['arms']:
            d=OUT/case/arm;r=load(d/'RECEIPT.json');row=dict(case_id=case,arm=arm,status=r['status'],final_NRMSE=None,raw_NRMSE=None,material_recovery_success='NO_ARBITRARY_GATE',
                original_Phi=r.get('original_Phi'),physical_source_RMS_relative_residual=r.get('physical_source_RMS_relative_residual'),
                CPU_seconds=r['phase2_incremental_CPU_seconds'],wall_seconds=r['phase2_incremental_wall_seconds'],counts=r.get('counts'),
                stop=r.get('optimizer'),dpcs=r.get('dpcs'),truth_access='OFFLINE_AFTER_ALL24_PUBLIC_SEAL')
            if r['status']=='SEALED_PUBLIC_PREDICTION':
                for file,key in [('raw.npz','raw_NRMSE'),('prediction.npz','final_NRMSE'),('checkpoint80_accepted.npz','checkpoint80_NRMSE')]:
                    if (d/file).exists():
                        with np.load(d/file,allow_pickle=False) as z:x=z['m'].copy()
                        row[key]=float(la.norm(x-truth)/norm)
                        if key=='final_NRMSE':final_material=x.copy()
                row['material_min']=float(final_material.min());row['material_max']=float(final_material.max())
                with np.load(d/'PAID_BUDGET_MATERIAL.npz',allow_pickle=False) as z:ms=z['m'].copy()
                costs=load(d/'PAID_BUDGET.json');offline=[]
                for j,(m,cost) in enumerate(zip(ms,costs)):offline.append(dict(index=j,chi_relative_L2_NRMSE=float(la.norm(m-truth)/norm),**cost))
                dump(OUT/'offline_curves'/f'{case}__{arm}.json',offline)
            scores.append(row)
    dump(OUT/'SCORES.json',dict(rows=scores,metric='||chi_hat-chi_reference||L2/||chi_reference||L2; NOT paperEq17 or scalaramplitudeRMSE',selection='paid originalPhiincumbents only; no GTselection',failed_arms_retained=True,score_new_LU=0,score_new_RHS=0))

def report(cfg):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figdir=ROOT/'FIGURES/vector_extension';figdir.mkdir(exist_ok=True)
    data=load(OUT/'SCORES.json')['rows'];cases=[x['case_id'] for x in cfg['cases']]
    fig,axes=plt.subplots(3,9,figsize=(20,7))
    for i,case in enumerate(cases):
        with np.load(TRUTH/case/'truth.npz') as z:truth=z['chi'].reshape((8,)*3)
        arrays=[truth]
        for arm in cfg['arms']:
            p=OUT/case/arm/'prediction.npz';arrays.append(np.load(p)['m'].reshape((8,)*3) if p.exists() else None)
        for j,x in enumerate(arrays):
            ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
            if x is None:ax.text(.5,.5,'FAILED',ha='center');continue
            ax.imshow(x[:,:,4].T,origin='lower',vmin=0,vmax=8,cmap='magma')
            if i==0:ax.set_title((['truth']+cfg['arms'])[j].replace('_','\n'),fontsize=8)
            if j==0:ax.set_ylabel(case,fontsize=9)
    fig.suptitle('Strict real 3D, all512unknowns: central-z slices (chi=epsilon-1)\nScaled/developmentDDA only;8->10crossgrid; allarms commonbox0..8');fig.tight_layout(rect=[0,0,.95,.91]);cax=fig.add_axes([.965,.2,.008,.6]);fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0,8),cmap='magma'),cax=cax,label='chi');fig.savefig(figdir/'material_slices.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(15,4.5))
    for ax,case in zip(axes,cases):
        for arm in cfg['arms']:
            p=OUT/'offline_curves'/f'{case}__{arm}.json'
            if not p.exists():continue
            rows=load(p);ax.plot([v['CPU_seconds'] for v in rows],[v['chi_relative_L2_NRMSE'] for v in rows],label=arm)
        ax.set_title(case);ax.set_xlabel('onlineCPU s');ax.set_ylabel('chi relative L2 NRMSE');ax.grid(alpha=.2)
    axes[-1].legend(fontsize=6);fig.suptitle('Offline-scored paid original-Phi incumbents; no GTiterate selection');fig.tight_layout();fig.savefig(figdir/'paid_cost_curves.png',dpi=160);plt.close(fig)
    lines=['# R9 Phase2：严格实材料三维矢量扩展','', '**本轮是scaled finite-DDA开发实验，不是原2019的3λ、30³/60³物理复现。** R8与PhaseI保持不变；所有24臂预测先封存，随后独立离线读GT评分。','', '三例均8³=512个独立实各向同性未知χ，完整1536维矢量场，60电偶极/60三分量接收，共10800复观测；所有方法共同χ盒[0,8]、初值.25和相同离散smooth/l2 prior。','', 'CIE单beta sourcewise采用6→1；Dual/global普通CIE current采用6→3。A/A+只替换材料读出，current保持与GlobalDual完全相同。原3D γ以频率/component键跨MF1→2精确transport，240固定PR、材料max12/A+两轮八色arc、共同300物理corrector/direct300，checkpoint80仅复制已付费incumbent。paper delta使用sqrt(mean ratio)作诊断，不触发提前停止，也不冒称paper Eq17 Err。','', 'fixed DPCS采用历史.1 threshold。新增qualified规则：旧10%固定能量case用(.1/.9)||Y||的确定性noise bound；新两例raw iid CN3%公开sigma并使用完整噪声球。两者均不覆盖fine/coarse离散模型误差，不是材料或连续物理证书。','', '|Case|Arm|Raw NRMSE|Final NRMSE|Checkpoint80 NRMSE|Residual|Wholeworker CPU s|LU|','|---|---|---:|---:|---:|---:|---:|---:|']
    for r in data:
        fmt=lambda v:'NOT_AVAILABLE' if v is None else f'{v:.5f}'
        lines.append(f"|{r['case_id']}|{r['arm']}|{fmt(r.get('raw_NRMSE'))}|{fmt(r.get('final_NRMSE'))}|{fmt(r.get('checkpoint80_NRMSE'))}|{fmt(r.get('physical_source_RMS_relative_residual'))}|{r['CPU_seconds']:.2f}|{(r.get('counts') or {}).get('physical_LU',0)}|")
    lines+=['','NRMSE为χ相对L2，与PhaseI幅度归一RMSE/成功门槛及paper Eq17不同；不跨池比较原始Phi。执行完成、optimizer停止、材料恢复和物理资格分别记录。','', '![实际材料](FIGURES/vector_extension/material_slices.png)','', '![已付费同成本轨迹](FIGURES/vector_extension/paid_cost_curves.png)','', '## 生成与物理资格','']
    for c in cases:
        rec=load(OUT/'generation'/f'{c}.json');lines.append(f"- {c}: fine10/coarse8 clean-data discrepancy {100*rec['crossgrid_clean_relative_difference']:.3f}%。")
    qual=load(OUT/'validation/RECEIPT.json');lines+=['', '独立Mie是εr9/χ8 sphere、.4λ箱、r=.145λ，固定24³子体素MG混合，8→10。不按结果改geometry/scale/预算，资格失败不会隐藏或停止其他反演。']
    for r in qual['rows']:lines.append(f"- grid{r['nside']}³, h={r['h']:.4f}λ, relative far-vector error={100*r['relative_Mie_far_vector_error']:.4f}%, linear residual={r['linear_relative_residual']:.3g}.")
    lines+=['', 'Mie error变化只是已测有限几何/contrast资格，不自动称连续精度PASS；MG不是精确cut-cell。其他宽度.6λ对象仍需按各自cross-grid差异解释，不能把.4λ sphere百分比移植过去。','', '## 成本与限制','', '所有新validation/generation/worker失败、setup、SVD、缓存命中、LU/source/adjointRHS和2进程队列费用在RESULTS/vector_extension与delegated/3d_extension/phase2_queue中保存。旧四块Y已复制到R9 public，clean/reference只在sealed_truth；其fine10模拟不重复收费。新两例各独立组装fine10和coarse8；sigma由生成器声明，solver只读publicY/采集几何/噪声metadata。','', 'χ8材料只用于新twocubes，cube side.12λ、center±.11λ、gap.10λ=两inverse cells；hollow沿用R8真正立方壳+杆，ε2.5。原四块ε2/3/4/5为.6λ缩小电尺寸改编。未知空间不缩成目标参数，不读取support。','', '本轮不训练NN、不安装后端、不运行原3λ实验。Direct300与联合240+300不是相同计算量；按真实paidCPU/LU轨迹比较，不能把上限退出称收敛或唯一恢复。']
    (ROOT/'VECTOR3D_EXTENSION.md').write_text('\n'.join(lines)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['all','prepare','worker','score','report']);p.add_argument('--case');p.add_argument('--arm');args=p.parse_args();cfg=load(CONFIG)
    if args.command=='worker':return worker(cfg,args.case,args.arm)
    if args.command=='prepare':prepare(cfg);return 0
    if args.command=='score':seal_and_score(cfg);return 0
    if args.command=='report':report(cfg);return 0
    if (OUT/'PHASE2_COMPLETE.json').exists():raise FileExistsError('CompletedPhase2preserved')
    started=utc();prepare(cfg);jobs=run_queue(cfg);seal_and_score(cfg);report(cfg)
    children=[]
    for case in [s['case_id'] for s in cfg['cases']]:
        for arm in cfg['arms']:children.append(load(OUT/case/arm/'RECEIPT.json'))
    dump(OUT/'PHASE2_COMPLETE.json',dict(status='COMPLETED' if all(j['returncode']==0 for j in jobs) else 'COMPLETED_WITH_FAILURES',started_utc=started,finished_utc=utc(),arms=24,
        parent_CPU_seconds=time.process_time()-ENTRY_CPU,children_CPU_seconds=sum(x['phase2_incremental_CPU_seconds'] for x in children),
        total_incremental_CPU_seconds=time.process_time()-ENTRY_CPU+sum(x['phase2_incremental_CPU_seconds'] for x in children),
        whole_phase_wall_seconds=time.perf_counter()-ENTRY_WALL,max_numeric_processes=2,threads_per_process=1,process_peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='parentCPUincludesimports/validation/generation/queue/score/figures; childCPUincludesallworkerimports/setup/solver/export; no nestedcostdoublecount; subprocesswall includesstartup',
        phaseI_modified=False,legacy_forward_reruns=0,truth_scoring='AFTER_ALL24_PUBLIC_SEAL',failures=[j for j in jobs if j['returncode']!=0]))
    return 0
if __name__=='__main__':
    with threadpool_limits(limits=1):raise SystemExit(main())
