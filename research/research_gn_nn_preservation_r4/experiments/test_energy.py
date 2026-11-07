"""Explicit post-seal TEST absolute-energy scoring. Never run during training.

All seal inventories, finite run plans and test observation identities are
validated before the first truth file is opened. No fitting/selection occurs.
"""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import argparse,csv,json,sys,time
from pathlib import Path
from collections import defaultdict
import numpy as np
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parents[1]
REPO=HERE.parents[1]
sys.path.insert(0,str(REPO))
from research.research_gn_nn_preservation_r4.src.energy import FixedEnergyGeometry,squared_norm
from research.research_z_alg_r4_diversity.src.physics import Model
from research.research_gn_nn_completion_r3.src.score import validate_stamp
from research.research_gn_nn_preservation_r4.experiments.e0_energy import writecsv

CORE=('total','P_in_strong19','P_in_weak557','S_strong9','S_weak39','S_visible48','S_null528','D509','B67','G_field','S_field')

def load_npz(path):
    with np.load(path,allow_pickle=False)as z:return {k:z[k]for k in z.files}

def validate_plan(seal_path,runs_path):
    """Public inventories only. No truth path is opened in this function."""
    seal=json.loads(Path(seal_path).read_text())
    if seal.get('status')!='ALL_PLANNED_PREDICTIONS_SEALED':raise ValueError('Complete prediction closure required')
    if Path(seal['study']).resolve()!=HERE.resolve():raise ValueError('This runner only scores the new R4 study')
    if set(seal.get('splits',[]))!={'test'}:raise ValueError('Explicit TEST-only seal required')
    if seal.get('max_observations')is not None:raise ValueError('Full TEST observation plan required; no test subset selection')
    plans=json.loads(Path(runs_path).read_text())
    if isinstance(plans,dict):plans=plans.get('plans',plans.get('runs'))
    if not isinstance(plans,list)or not plans:raise ValueError('Explicit nonempty frozen run list required')
    plans=[dict(run=p,arm=p,seed=None,passes=None)if isinstance(p,str)else dict(p)for p in plans]
    runs=[p['run']for p in plans]
    if len(set(runs))!=len(runs)or set(runs)!=set(seal['runs']):raise ValueError('Explicit run list must exactly equal the seal run plan')
    validate_stamp(seal['manifest'])
    manifest=json.loads(Path(seal['manifest']['path']).read_text())
    records=manifest['records']
    if any(r['split']!='test'for r in records):raise ValueError('Sealed manifest contains non-TEST records')
    online=json.loads((HERE/'data/manifests/online.json').read_text())['records']
    expected={r['obs_id']for r in online if r['split']=='test'}
    ids=[r['obs_id']for r in records]
    if len(ids)!=len(set(ids))or set(ids)!=expected or set(ids)!=set(seal['observation_ids']):raise ValueError('TEST observation identities do not match the frozen full public plan')
    items={}
    for item in seal['receipts']:
        key=(item['run'],item['obs_id'])
        if key in items:raise ValueError('Duplicate sealed run-observation receipt')
        items[key]=item
        if item['disposition']=='SEALED':
            validate_stamp(item['json']);validate_stamp(item['npz'])
        elif item['disposition']!='NOT_RUN':raise ValueError('Unknown sealed receipt disposition')
    if set(items)!={(run,obs)for run in runs for obs in ids}:raise ValueError('Incomplete finite closure inventory')
    return seal,plans,records,items

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--seal',type=Path,required=True)
    p.add_argument('--runs-json',type=Path,required=True,help='Explicit frozen list, or {plans:[{run,arm,seed,passes}]}')
    p.add_argument('--output',type=Path,default=HERE/'results/test_energy')
    a=p.parse_args();start=time.perf_counter()
    seal,plans,records,items=validate_plan(a.seal,a.runs_json)
    # First truth access can occur only below this completed validation boundary.
    out=a.output;out.mkdir(parents=True,exist_ok=True)
    ops=load_npz(HERE/'data/operators.npz')
    geom=FixedEnergyGeometry.from_cache(REPO/'research/research_gn_green_error_20261106/results/fixed_spectra.npz',G=ops['G'],S=ops['S'],sensitivity=True)
    model=Model(ops);truth_cache={};notrun=[];residuals=[];realization=[];stats=defaultdict(int);defects=[]
    aggregate=defaultdict(lambda:dict(n=0,sums=defaultdict(float)))
    boot=defaultdict(lambda:defaultdict(list))
    fp=(out/'per_case_energy.csv').open('w',newline='');writer=None
    meta_keys=('run','arm','seed','passes','obs_id','parent_id','group','family','noise','amplitude')
    group_keys=('run','arm','seed','passes','group','family','noise','amplitude','stage','quantity','projection')
    bootstrap_keys=('run','arm','seed','passes','noise','cohort','stage','quantity','projection')
    def emit(meta,quantity,stage,before,after):
        nonlocal writer
        row=dict(meta,quantity=quantity,stage=stage)
        for band,values in geom.transition(before,after).items():
            for name,value in values.items():row[f'{band}__{name}']=value
            key=tuple(row[k]for k in group_keys[:-1])+(band,)
            group=aggregate[key];group['n']+=1
            for name,value in values.items():group['sums'][name]+=value
            if band in CORE:
                for cohort in ('ALL',meta['group']):
                    br=dict(row,cohort=cohort,projection=band)
                    bk=tuple(br[k]for k in bootstrap_keys)
                    boot[bk][meta['parent_id']].append(values['delta'])
        if writer is None:writer=csv.DictWriter(fp,fieldnames=list(row));writer.writeheader()
        writer.writerow(row)
    with threadpool_limits(limits=1):
        for oi,r in enumerate(sorted(records,key=lambda z:z['obs_id'])):
            obs=r['obs_id'];public=load_npz(HERE/'data'/r['public_path']);Y=public['Y']
            active=[p for p in plans if items[p['run'],obs]['disposition']=='SEALED']
            if not active:
                notrun.extend(dict(run=p['run'],obs_id=obs,reason='SEALED_NOT_RUN')for p in plans);continue
            truth_key=r['truth_path']
            if truth_key not in truth_cache:
                truth=load_npz(HERE/'data'/truth_key) # explicit first offline TEST truth boundary
                mt=np.asarray(truth['m']);st=model.forward(mt,label='TEST_energy_unique_GT')
                truth_cache[truth_key]=(mt,st['j'].copy(),st['y'].copy());defects.append(st['state_relative_residual'])
                stats['unique_GT_solves']+=1
            mt,jt,yt=truth_cache[truth_key];eta=Y-yt;n2=squared_norm(eta);y2=squared_norm(Y)
            currents={};anchor_done=False
            def physical(m,label):
                key=np.asarray(m,dtype=float).tobytes()
                if key in currents:stats['material_cache_hits']+=1;return currents[key]
                st=model.forward(np.asarray(m,float),label=label);defects.append(st['state_relative_residual'])
                currents[key]=(st['j'].copy(),st['y'].copy());stats['unique_material_solves']+=1
                return currents[key]
            for plan in plans:
                run=plan['run'];item=items[run,obs]
                if item['disposition']!='SEALED':notrun.append(dict(run=run,obs_id=obs,reason='SEALED_NOT_RUN'));continue
                arrays=load_npz(item['npz']['path']);receipt=json.loads(Path(item['json']['path']).read_text())
                ma=np.asarray(arrays['anchor']);ja,ya=physical(ma,'TEST_energy_anchor')
                meta={k:r.get(k)for k in ('obs_id','parent_id','group','family','noise','amplitude')}
                meta.update(run=run,arm=plan.get('arm',run),seed=plan.get('seed'),passes=plan.get('passes'))
                raw_receipt=next((c for c in receipt.get('candidates',[])if c['label']=='raw_1'),None)
                raw=None
                if raw_receipt is not None:
                    raw=np.asarray(arrays['candidate_materials'][raw_receipt['index']])
                    jr,yr=physical(raw,'TEST_energy_frozen_verified_raw_1')
                    emit(meta,'material_complex_Green_embedding','anchor_to_raw_1',ma-mt,raw-mt)
                    emit(meta,'Maxwell_current','anchor_to_raw_1',ja-jt,jr-jt)
                else:notrun.append(dict(run=run,obs_id=obs,reason='NO_FROZEN_RAW_1_CANDIDATE_raw_energy_NOT_RUN'))
                if 'j_aux'in arrays:
                    emit(meta,'auxiliary_current_head_OFFLINE_not_Maxwell_solution','anchor_to_auxiliary_head',ja-jt,np.asarray(arrays['j_aux'])-jt)
                    if raw is not None:
                        jaux=np.asarray(arrays['j_aux']);qaux=ops['E0']+ops['G']@jaux
                        state_defect=jaux-raw[:,None]*qaux
                        gap=jr-jaux;aux_receiver=ops['S']@(jaux-ja);receiver_gap=ops['S']@gap
                        realization.append(dict(meta,state_defect_energy=squared_norm(state_defect),
                            current_realization_gap_energy=squared_norm(gap),
                            internal_realization_gap_energy=squared_norm(ops['G']@gap),
                            aux_receiver_update_energy=squared_norm(aux_receiver),
                            receiver_realization_gap_energy=squared_norm(receiver_gap),
                            signed_receiver_cross=float(2*np.vdot(aux_receiver,receiver_gap).real),
                            physical_receiver_update_energy=squared_norm(yr-ya),
                            definition='F(raw)-F(anchor)=S(jaux-janchor)+S(jphysical_raw-jaux); current realization gap=-L_raw^-1 state_defect',
                            scope='NN matched auxiliary/readout when applicable; baseline raw incumbents need their own current pairing interpretation'))
                final=np.asarray(arrays['top1'])if'top1'in arrays else None
                if final is not None:
                    jf,yf=physical(final,'TEST_energy_frozen_top1')
                    for qty,be,fe in [('material_complex_Green_embedding',ma-mt,final-mt),('Maxwell_current',ja-jt,jf-jt)]:emit(meta,qty,'anchor_to_final_top1',be,fe)
                    if raw is not None:
                        emit(meta,'material_complex_Green_embedding','raw_1_to_final_pipeline',raw-mt,final-mt)
                        emit(meta,'Maxwell_current','raw_1_to_final_pipeline',jr-jt,jf-jt)
                    stats['final_top1_available']+=1
                else:notrun.append(dict(run=run,obs_id=obs,reason='NO_FROZEN_FEASIBLE_TOP1_final_energy_NOT_RUN'))
                states=[('anchor',ya)]
                if raw is not None:states.append(('raw_1',yr))
                if final is not None:states.append(('final_top1',yf))
                for stage,yp in states:
                    resid=yp-Y
                    residuals.append(dict(meta,stage=stage,true_S_field_error_energy=squared_norm(yp-yt),observed_misfit_energy=squared_norm(resid),noise_energy=n2,signed_residual_noise_cross=float(2*np.vdot(resid,eta).real),observed_data_relative_residual=float(np.sqrt(squared_norm(resid)/y2)),identity='Fpred-FGT = (Fpred-Y)+eta; eta=Y-FGT'))
                stats['sealed_model_observations']+=1
            if oi%8==0 or oi==len(records)-1:
                fp.flush();progress=dict(observations_completed=oi+1,total_observations=len(records),counts=dict(stats),wall_seconds=time.perf_counter()-start)
                (out/'PROGRESS.json').write_text(json.dumps(progress,indent=2));print(json.dumps(progress),flush=True)
    fp.close();rows=[]
    for key,v in aggregate.items():
        row=dict(zip(group_keys,key));row['cases']=v['n']
        for name,value in v['sums'].items():row[name+'_sum']=value;row[name+'_mean']=value/v['n']
        rows.append(row)
    writecsv(out/'grouped_absolute_energy.csv',rows)
    rng=np.random.default_rng(2026110705);brows=[]
    for key,mothers in boot.items():
        values=np.array([np.mean(v)for _,v in sorted(mothers.items())]);n=len(values)
        means=values[rng.integers(0,n,(1000,n))].mean(1)
        brows.append(dict(zip(bootstrap_keys,key),mothers=n,delta_mean=float(values.mean()),delta_ci_low=float(np.quantile(means,.025)),delta_ci_high=float(np.quantile(means,.975)),mother_improvement_fraction=float(np.mean(values<0)),bootstrap_repetitions=1000))
    writecsv(out/'mother_paired_bootstrap.csv',brows);writecsv(out/'data_residual_relation.csv',residuals);writecsv(out/'NOT_RUN.csv',notrun)
    writecsv(out/'current_realization_gap.csv',realization)
    receipt=dict(status='COMPLETE_FINITE_SEALED_PLAN_AUDITED',seal_path=str(a.seal.resolve()),runs_json=str(a.runs_json.resolve()),
        planned_models=len(plans),planned_observations=len(records),counts=dict(stats),physical_costcounts=dict(model.counts),
        max_full_VIE_equation_defect=max(defects,default=None),wall_seconds=time.perf_counter()-start,
        BLAS_threads=1,CPU_workers=1,truth_access='NEW_TEST_ONLY_AFTER_FULL_INVENTORY_VALIDATION',
        raw_definition='actual frozen raw_1 candidate material, including its original box clipping',
        auxiliary_scope='separately named head discrepancies, never substituted for nonlinear Maxwell currents',
        material_scope='complex Green coordinate embedding only, not real material observability',
        no_new_optimization=True,no_new_training=True,geometry=geom.metadata,
        scientific_scope='absolute errors and mother-paired changes; finite-plan NOT_RUN entries retained; no test-based method choice')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
