"""Frozen CM nine-arm finite job, two public Y and one resident G/E/S/Vr.

No material truth is loaded. Matched Y reuses already paid clean observations;
all nine terminal predictions must be sealed before separate offline scoring.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import json,sys,time,traceback,gc,platform
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np
import scipy
from scipy import linalg as la
import psutil,torch
from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel
from research.research_max_struct_r9.src.vector_original_cuda import CudaOriginalModel,CudaReceiverGeometry,CudaReceiverSubspace,cuda_auxiliary_solve
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d'

class ControlledStop(SystemExit):pass
def dump(p,v):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def counts_delta(model,before):return {k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}

def public_data(protocol):
    wall,cput=time.perf_counter(),time.process_time()
    old=np.load(OUT/'public/observations.npz');truth_fields=np.load(OUT/'sealed_truth/truth.npz')
    # Read only already paid clean/noise observation fields, never material GT.
    fine=truth_fields['Y_clean'];noise=truth_fields['noise'];coarse=np.load(OUT/'sealed_truth/coarse_truth_response.npz')['prediction']
    scale=float(la.norm(coarse)/la.norm(fine));matched_noise=noise*scale
    Ys={'cross_grid':old['Y'].copy(),'matched_grid':coarse+matched_noise}
    sigmas={'cross_grid':float(old['sigma']),'matched_grid':float(old['sigma'])*scale}
    clean={'cross_grid':fine,'matched_grid':coarse};noises={'cross_grid':noise,'matched_grid':matched_noise};result={}
    for cohort in protocol['cohorts']:
        name=cohort['name'];dest=OUT/'cohorts'/name;file=dest/'public/observations.npz'
        if file.exists():raise RuntimeError('Refusing overwrite of frozen cohort '+name)
        file.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(file,Y=Ys[name],sigma=sigmas[name],noise_radius=sigmas[name]*np.sqrt(Ys[name].size),
            data_id=cohort['data_id'],inverse_side=30,width_lambda=3.,material_box=np.array([0.,4.]))
        private=dest/'sealed_observation_reference';private.mkdir(exist_ok=True)
        np.savez_compressed(private/'fields.npz',Y_clean=clean[name],noise=noises[name])
        stats=dict(data_id=cohort['data_id'],scope=protocol['scope'],sigma=sigmas[name],
            noise_relative_realized=float(la.norm(noises[name])/la.norm(clean[name])),noise_relative_expected=.1,
            seed=2026110619,shared_noise_direction=True,noise_rule=protocol['matched_noise'],
            Y_shape=list(Ys[name].shape),new_forward_generation_RHS=0,material_GT_in_public_inputs=False,
            clean_source='Previously paid CM60 clean prediction' if name=='cross_grid' else 'Previously paid CM30 clean prediction',
            cross_Y_unchanged=bool(np.array_equal(Ys[name],old['Y'])) if name=='cross_grid' else None)
        dump(dest/'public/manifest.json',stats);result[name]=dict(Y=Ys[name],sigma=sigmas[name],stats=stats)
    dump(OUT/'NINE_DATA_PREPARATION_COST.json',dict(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,
        forward_generation_RHS=0,noise_scale_coarse_over_fine=scale,
        shared_direction_relative_difference=float(la.norm(noise/la.norm(noise)-matched_noise/la.norm(matched_noise))),
        cohorts={n:v['stats'] for n,v in result.items()},additive=False,included_in='run_nine_CM_PROCESS_COST.json'))
    return result

def main():
    psutil.Process().cpu_affinity([12,14]);torch.set_num_threads(1);torch.set_num_interop_threads(1)
    config=json.loads((BASE/'configs/original_vector3d.json').read_text());protocol=json.loads((BASE/'configs/original_cm_nine.json').read_text())
    if (OUT/'SEALED_ALL_NINE.json').exists():raise RuntimeError('Nine endpoints already sealed')
    dump(OUT/'ACTIVE_CUDA_JOB.json',dict(pid=os.getpid(),mode='CM_nine',process_started_epoch=psutil.Process().create_time(),affinity=[12,14]))
    model=None;status='FAILED';failure=None;failures=[];cache={};cache_rows=[];terminals=[]
    try:
        if (OUT/'STOP_CURRENT_JOB_NINE').exists():raise ControlledStop('Explicit stop marker')
        datasets=public_data(protocol)
        # Frozen CM constructor only: no IGT import or substitution.
        model=CudaOriginalModel(OriginalVectorModel(30),fft_layout='packed',block_rhs=16)
        geometry=CudaReceiverGeometry(model,.1)
        shared=dict(CPU_seconds=time.process_time(),wall_seconds=time.time()-psutil.Process().create_time(),
            setup=model.setup,CPU_operator_setup_counts=dict(model.cpu.counts),counts=dict(model.counts),receiver_geometry=geometry.info,
            source='One CM G/E/S and SGram/EVD/Vr for both Y/all9arms',additive=False,included_in='run_nine_CM_PROCESS_COST.json',
            includes_startup_and_data_preparation=True,amortization_arms=9)
        dump(OUT/'NINE_SHARED_OPERATOR_COST.json',shared)
        for cohort in protocol['cohorts']:
            name=cohort['name'];before=dict(model.counts);wall,cput=time.perf_counter(),time.process_time()
            Y=model.array(datasets[name]['Y']);receiver=CudaReceiverSubspace(model,Y,.1,geometry=geometry)
            setup=dict(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,counts=counts_delta(model,before),
                receiver=receiver.info,data_id=cohort['data_id'],additive=False,included_in='run_nine_CM_PROCESS_COST.json',
                shared_geometry_not_recomputed=True)
            dump(OUT/'cohorts'/name/'RECEIVER_Y_SETUP_COST.json',setup)
            datasets[name].update(Y_cuda=Y,receiver=receiver,setup=setup)

        def physical_cached(m,Yhost,label,progress):
            mhost=m.cpu().numpy();key=mhost.tobytes();model.counts['physical_value_calls']+=1
            if key in cache:
                item=cache[key];model.counts['physical_prediction_cache_hits']+=1;hit=True;actual_rhs=0
            else:
                wall,cput=time.perf_counter(),time.process_time();before=dict(model.counts);first=len(model.solve_records)
                if np.all(mhost==0):
                    prediction=np.zeros((180,60),complex);model.counts['analytic_air_zero_predictions']+=1
                else:prediction=model.physical(m,label=label,progress=progress)
                actual_rhs=len(model.solve_records)-first
                item=dict(id='physical%03d'%len(cache),prediction=prediction,first_label=label,uses=[],actual_first_RHS=actual_rhs,
                    evaluation_CPU_seconds=time.process_time()-cput,evaluation_wall_seconds=time.perf_counter()-wall,
                    counts=counts_delta(model,before),max_relative_linear_residual=max([r.get('relative_residual',0.) for r in model.solve_records[first:]] or[0.]))
                cache[key]=item;cache_rows.append(item);hit=False
            item['uses'].append(label)
            dump(OUT/'NINE_PHYSICAL_CACHE.json',dict(scope='Exact material bytes, same fixed G/E/S; no Y or GT criterion',
                entries=[{k:v for k,v in item.items() if k!='prediction'} for item in cache_rows],
                additive=False,included_in='run_nine_CM_PROCESS_COST.json'))
            pred=item['prediction'];loss=float(np.mean(np.sum(abs(pred-Yhost)**2,axis=0)/np.sum(abs(Yhost)**2,axis=0)))
            return dict(prediction=pred.copy(),data_loss_mean=loss,relative_source_RMS=float(np.sqrt(loss)),
                physical_cache_id=item['id'],physical_cache_hit=hit,actual_physical_RHS=actual_rhs,
                reused_prediction_first_label=item['first_label'],max_relative_linear_residual=item['max_relative_linear_residual'])

        for cohort in protocol['cohorts']:
            name=cohort['name'];data=datasets[name];narm=len(cohort['arms'])
            for arm in cohort['arms']:
                dest=OUT/'cohorts'/name/'inversions'/arm
                if (dest/'SEALED.json').exists():raise RuntimeError('Duplicate frozen endpoint')
                dest.mkdir(parents=True,exist_ok=True);before=dict(model.counts);start_records=len(model.solve_records)
                awall,acpu=time.perf_counter(),time.process_time();endpoints=[];armstatus='FAILED';error=None
                overheadCPU=shared['CPU_seconds']/9+data['setup']['CPU_seconds']/narm
                overheadWall=shared['wall_seconds']/9+data['setup']['wall_seconds']/narm
                dump(dest/'FROZEN_CONFIG.json',dict(config,cohort=name,data_id=cohort['data_id'],arm=arm,
                    operator='CM+RR full30^3',truth_access=False,primary_endpoint='stage2',protocol='CM_NINE_PROTOCOL_AMENDMENT_ZH.md'))
                try:
                    def callback(r,m):
                        r.update(cumulative_CPU_seconds=time.process_time()-acpu+overheadCPU,
                            cumulative_wall_seconds=time.perf_counter()-awall+overheadWall,
                            process_CPU_seconds=time.process_time(),runtime_counts=dict(model.counts))
                        with (dest/'history.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                        if r['iteration']==1 or r['iteration']%25==0 or (r['coefficient_delta_RMS'] is not None and r['coefficient_delta_RMS']<.001):
                            print(json.dumps(dict(cohort=name,arm=arm,**r)),flush=True)
                        if (OUT/'STOP_CURRENT_JOB_NINE').exists():
                            np.savez_compressed(dest/'controlled_stop_material.npz',material=m.cpu().numpy())
                            raise ControlledStop('Explicit stop at fully accounted iteration boundary')
                    def endpoint(epname,m,stage):
                        np.savez_compressed(dest/(epname+'.npz'),material=m.cpu().numpy())
                        label=name+'_'+arm+'_'+epname
                        def progress(r):
                            with (dest/'physical_progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                            if r['source']%20==0:print(json.dumps(dict(cohort=name,arm=arm,**r)),flush=True)
                        value=physical_cached(m,data['Y'],label,progress)
                        np.savez_compressed(dest/(epname+'_prediction.npz'),prediction=value.pop('prediction'))
                        endpoints.append(dict(name=epname,stage=stage,**value,
                            cumulative_CPU_seconds=time.process_time()-acpu+overheadCPU,
                            cumulative_wall_seconds=time.perf_counter()-awall+overheadWall))
                        dump(dest/'ENDPOINTS.json',endpoints)
                    result=cuda_auxiliary_solve(model,data['receiver'],data['Y_cuda'],arm,config,callback,endpoint)
                    monotonic=[]
                    if arm in ('sourceCIE','LSIEsource'):
                        for stage in range(2):
                            values=np.array([h['auxiliary_sum'] for h in result['history'] if h['stage']==stage]);d=np.diff(values)
                            monotonic.append(dict(stage=stage,max_absolute_increase=float(max(0,d.max())) if len(d) else 0.,
                                maximum_signed_change=float(d.max()) if len(d) else 0.,max_relative_increase=float(max(0,np.max(d/abs(values[:-1])))) if len(d) else 0.))
                    summary=dict(cohort=name,arm=arm,data_id=cohort['data_id'],stages=result['stages'],final=endpoints[-1],
                        receiver_rank=geometry.rank,within_stage_monotonicity_diagnostic=monotonic,prior=0,rho=1,
                        primary_endpoint='Final stage, no GT choice',physics_cache_scope='Actual shared queue cost; not standalone algorithm speed')
                    dump(dest/'SUMMARY.json',summary);dump(dest/'SEALED.json',dict(status='PREDICTIONS_SEALED_BEFORE_ANY_NINE_GT_SCORE',**summary));armstatus='COMPLETE'
                    terminals.append(dict(cohort=name,arm=arm,status=armstatus,seal=str(dest.relative_to(OUT)/'SEALED.json')))
                except BaseException as exc:
                    error=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());failures.append(dict(cohort=name,arm=arm,**error))
                    if isinstance(exc,(KeyboardInterrupt,SystemExit)):raise
                    dump(dest/'TERMINAL_FAILURE.json',dict(status='TERMINAL_NUMERICAL_FAILURE_NO_GT_SELECTION',failure=error))
                finally:
                    torch.cuda.synchronize()
                    dump(dest/'COST.json',dict(status=armstatus,platform='XINAN',additive=False,included_in='run_nine_CM_PROCESS_COST.json',
                        CPU_seconds=time.process_time()-acpu,wall_seconds=time.perf_counter()-awall,
                        setup_amortized_CPU_seconds=time.process_time()-acpu+overheadCPU,
                        setup_amortized_wall_seconds=time.perf_counter()-awall+overheadWall,
                        shared_operator_CPU_share=shared['CPU_seconds']/9,Y_setup_CPU_share=data['setup']['CPU_seconds']/narm,
                        counts=counts_delta(model,before),solve_records=model.solve_records[start_records:],
                        scope='Actual shared queue with exact-m prediction cache; standalone runtime not measured',failure=error))
                gc.collect()
        status='COMPLETE_WITH_ARM_FAILURES' if failures else 'COMPLETE'
        if len(terminals)==9:
            dump(OUT/'SEALED_ALL_NINE.json',dict(status='ALL_NINE_FINAL_ENDPOINTS_SEALED_BEFORE_GT_SCORING',terminals=terminals,
                protocol=protocol,truth_scoring_has_not_run=True))
    except BaseException as exc:
        status='CONTROLLED_STOP' if isinstance(exc,ControlledStop) else 'FAILED'
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());raise
    finally:
        if model:torch.cuda.synchronize()
        ct=psutil.Process().cpu_times()
        dump(OUT/'run_nine_CM_PROCESS_COST.json',dict(status=status,platform='XINAN',additive=True,included_in=None,
            CPU_seconds=ct.user+ct.system,wall_seconds=time.time()-psutil.Process().create_time(),
            scope='Whole process including imports/data prep/shared G/E/S/Vr/two Y fits/all9arms/endpoint validation/cache/I/O',
            counts=dict(model.counts) if model else{},CPU_operator_setup_counts=dict(model.cpu.counts) if model else{},
            solve_records=model.solve_records if model else[],terminals=terminals,arm_failures=failures,failure=failure,
            peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated() if model else 0,
            peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved() if model else 0,
            peak_RSS_bytes=psutil.Process().memory_info().peak_wset,process_affinity=psutil.Process().cpu_affinity(),
            platform_version=platform.platform(),numpy_version=np.__version__,scipy_version=scipy.__version__,torch_version=torch.__version__))

if __name__=='__main__':main()
