"""Offline E0: absolute historical GN/NN error energies; no new training."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import argparse,csv,json,sys,time
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parents[1]
REPO=HERE.parents[1]
sys.path.insert(0,str(REPO))
from research.research_gn_nn_preservation_r4.src.energy import FixedEnergyGeometry,squared_norm
from research.research_z_alg_r4_diversity.src.physics import Model

META=('obs_id','parent_id','family','distribution','amplitude','level','noise','run','arm','seed','passes','stage','quantity')
CORE=('total','P_in_strong19','P_in_weak557','S_strong9','S_weak39','S_visible48','S_null528','D509','B67','G_field','S_field')

def writecsv(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--max-minutes',type=float,default=35.)
    a=p.parse_args(); start=time.perf_counter(); deadline=start+60*a.max_minutes
    out=HERE/'results/e0';out.mkdir(parents=True,exist_ok=True)
    figs=HERE/'figures/e0';figs.mkdir(parents=True,exist_ok=True)
    evidence=HERE/'delegated/energy';evidence.mkdir(parents=True,exist_ok=True)
    old=REPO/'research/research_gn_nn_completion_r3'
    protocol=json.loads((HERE/'audit/PROTOCOL.json').read_text())
    seal=json.loads((old/'evaluation/CONFIRMATION_SEALED.json').read_text())
    if seal['status']!='ALL_PLANNED_PREDICTIONS_SEALED':raise RuntimeError('Historical predictions not sealed')
    plans=json.loads((old/'audit/CONFIRMATION_PLAN_FROZEN.json').read_text())['plans']
    records=json.loads((old/'evaluation/CONFIRMATION_SEALED_MANIFEST.json').read_text())['records']
    records=sorted(records,key=lambda r:(r['level'],r['parent_id'],r['noise']))
    public=dict(np.load(old/'data/operators.npz',allow_pickle=False))
    geom=FixedEnergyGeometry.from_cache(REPO/'research/research_gn_green_error_20261106/results/fixed_spectra.npz',G=public['G'],S=public['S'])
    (out/'GEOMETRY.json').write_text(json.dumps(geom.metadata,indent=2))
    model=Model(public);gt_cache={};counts=defaultdict(int);defects=[];notrun=[]
    # Selection is declared from public historical metadata, before errors read.
    chosen=defaultdict(set)
    for r in records:
        if r['level']=='low':
            key=(r['family'],r['distribution'])
            if len(chosen[key])<2:chosen[key].add(r['parent_id'])
    diag_ids={r['obs_id'] for r in records if r['level']=='low' and r['parent_id'] in chosen[r['family'],r['distribution']]}
    (out/'J_DIAGNOSTIC_SELECTION.json').write_text(json.dumps(dict(rule='lexicographic first two low-amplitude parents per family × ID/OOD, both noise observations',observations=sorted(diag_ids)),indent=2))
    summary=defaultdict(lambda:dict(n=0,sums=defaultdict(float)))
    # Mother-paired bootstrap applies directly to after-before differences.
    boot=defaultdict(lambda:defaultdict(list));jrows=[];case_rows=0
    fp=(out/'per_case_energy.csv').open('w',newline='');writer=None
    def emit(meta,quantity,stage,before,after):
        nonlocal writer,case_rows
        trans=geom.transition(before,after)
        row=dict(meta,quantity=quantity,stage=stage)
        for band,vals in trans.items():
            for metric,value in vals.items():row[f'{band}__{metric}']=value
            keys=('run','arm','seed','passes','noise','amplitude','level','stage','quantity')
            key=tuple(row[k] for k in keys)+(band,)
            q=summary[key];q['n']+=1
            for k,v in vals.items():q['sums'][k]+=v
            if band in CORE:
                bk=tuple(row[k] for k in ('run','arm','seed','passes','noise','level','stage','quantity'))+(band,)
                boot[bk][row['parent_id']].append(vals['delta'])
        if writer is None:writer=csv.DictWriter(fp,fieldnames=list(row));writer.writeheader()
        writer.writerow(row);case_rows+=1
    with threadpool_limits(limits=1):
        for oi,r in enumerate(records):
            if time.perf_counter()>deadline:
                notrun.extend(dict(obs_id=z['obs_id'],run=p['run'],reason='E0_LOCAL_TIME_LIMIT') for z in records[oi:] for p in plans)
                break
            obs=r['obs_id'];key=r['truth_path']
            if key not in gt_cache:
                with np.load(old/'data'/key,allow_pickle=False) as z:mt=np.array(z['m'])
                st=model.forward(mt,label='E0_unique_GT');gt_cache[key]=(mt,st['j'].copy())
                defects.append(st['state_relative_residual']);counts['unique_GT_solves']+=1
            mt,jt=gt_cache[key]
            with np.load(old/'anchors'/f'{obs}.npz',allow_pickle=False) as z:ma=np.array(z['m0'])
            ast=model.forward(ma,label='E0_anchor');ja=ast['j'].copy();defects.append(ast['state_relative_residual']);counts['anchor_solves']+=1
            JV=None
            if obs in diag_ids:
                J=model.jacobian(ast)
                _,s,vh=svd(J,full_matrices=False,check_finite=False)
                JV=vh[s>=s[0]*.1].T.copy();counts['real_anchor_J_diagnostics']+=1
            materials={ma.tobytes():ja}
            def physical(m):
                k=np.asarray(m,dtype=float).tobytes()
                if k in materials:counts['within_observation_material_cache_hits']+=1;return materials[k]
                st=model.forward(np.asarray(m,float),label='E0_frozen_material')
                materials[k]=st['j'].copy();defects.append(st['state_relative_residual']);counts['raw_or_final_unique_solves']+=1
                return materials[k]
            for plan in plans:
                path=old/'evaluation'/plan['run']/f'{obs}.npz'
                meta={k:r[k] for k in META[:7]};meta.update({k:plan[k] for k in ('run','arm','seed','passes')})
                if not path.exists():notrun.append(dict(obs_id=obs,run=plan['run'],reason='FROZEN_EVALUATION_MISSING'));continue
                with np.load(path,allow_pickle=False) as z:
                    raw=np.array(z['m_proposal_unclipped']); final=np.array(z['top1']) if 'top1' in z else None
                    aux=np.array(z['j_aux']) if 'j_aux' in z else None
                jr=physical(raw)
                emit(meta,'material_complex_Green_embedding','anchor_to_raw_unclipped',ma-mt,raw-mt)
                emit(meta,'Maxwell_current','anchor_to_raw_unclipped',ja-jt,jr-jt)
                counts['frozen_model_observations']+=1
                if aux is not None:
                    emit(meta,'auxiliary_current_head_OFFLINE_not_Maxwell_solution','anchor_to_auxiliary_head',ja-jt,aux-jt)
                    counts['auxiliary_heads']+=1
                if final is not None:
                    jf=physical(final)
                    for qty,be,re,fe in [('material_complex_Green_embedding',ma-mt,raw-mt,final-mt),('Maxwell_current',ja-jt,jr-jt,jf-jt)]:
                        emit(meta,qty,'anchor_to_final_top1',be,fe)
                        emit(meta,qty,'raw_to_final_correction',re,fe)
                    counts['final_top1_available']+=1
                else:
                    notrun.append(dict(obs_id=obs,run=plan['run'],reason='NO_FROZEN_PHYSICALLY_FEASIBLE_TOP1_final_energy_NOT_RUN'))
                if JV is not None:
                    for stage,m in [('anchor',ma),('raw_unclipped',raw),('final_top1',final)]:
                        if m is None:continue
                        e=m-mt;total=squared_norm(e);strong=squared_norm(JV.T@e)
                        jrows.append(dict(meta,stage=stage,real_J_strong_rank=JV.shape[1],material_error_energy=total,J_strong_error_energy=strong,J_weak_error_energy=total-strong,definition='real analytic material Jacobian at frozen GN anchor; threshold sigma/sigma_max >= .1'))
            if oi%8==0 or oi==len(records)-1:
                fp.flush();(out/'PROGRESS.json').write_text(json.dumps(dict(observations_completed=oi+1,total_observations=len(records),counts=dict(counts),elapsed_seconds=time.perf_counter()-start),indent=2))
                print(json.dumps(dict(observations=oi+1,rows=case_rows,seconds=round(time.perf_counter()-start,1),LU=model.counts['LU'])),flush=True)
    fp.close()
    groupkeys=('run','arm','seed','passes','noise','amplitude','level','stage','quantity','projection')
    srows=[]
    for key,v in sorted(summary.items(),key=lambda z:str(z[0])):
        row=dict(zip(groupkeys,key));row['cases']=v['n']
        for name,value in v['sums'].items():row[name+'_sum']=value;row[name+'_mean']=value/v['n']
        srows.append(row)
    writecsv(out/'grouped_absolute_energy.csv',srows)
    rng=np.random.default_rng(2026110704);brows=[]
    bkeys=('run','arm','seed','passes','noise','level','stage','quantity','projection')
    for key,byparent in sorted(boot.items(),key=lambda z:str(z[0])):
        values=np.array([np.mean(v) for _,v in sorted(byparent.items())]);n=len(values)
        draws=rng.integers(0,n,size=(1000,n));means=np.mean(values[draws],axis=1)
        brows.append(dict(zip(bkeys,key),mothers=n,delta_mean=float(np.mean(values)),delta_ci_low=float(np.quantile(means,.025)),delta_ci_high=float(np.quantile(means,.975)),mother_improvement_fraction=float(np.mean(values<0)),bootstrap_repetitions=1000,bootstrap_unit='mother parent, paired after-before absolute energy; seeds not pooled'))
    writecsv(out/'mother_paired_bootstrap.csv',brows);writecsv(out/'real_J_diagnostic.csv',jrows);writecsv(out/'NOT_RUN.csv',notrun)
    receipt=dict(status='COMPLETE' if counts['frozen_model_observations']==len(records)*len(plans) else 'PARTIAL',
        historical_predictions_sealed=True,new_training=False,new_test_truth_access=False,
        historical_exposed_truth_scope='R3 test only, offline frozen-prediction audit',
        requested_observations=len(records),requested_models=len(plans),expected_model_observations=len(records)*len(plans),
        low_parents=64,stress_parents=32,counts=dict(counts),physical_costcounts=dict(model.counts),
        max_full_VIE_equation_defect=max(defects,default=None),wall_seconds=time.perf_counter()-start,
        BLAS_threads=1,CPU_workers=1,raw_definition='frozen m_proposal_unclipped, no new clipping',
        final_definition='frozen physically selected top1 only; absent top1 explicitly NOT_RUN',
        auxiliary_definition='auxiliary-head discrepancy to GT current, never called true Maxwell current error',
        bootstrap='paired absolute energy differences with 1000 mother resamples; seed-specific; low/stress separate',
        projection_cache=str(REPO/'research/research_gn_green_error_20261106/results/fixed_spectra.npz'),
        physics_source='unchanged research/research_z_alg_r4_diversity/src/physics.py Model')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    residual_relation(old,out,public['S'])
    plots_and_findings(srows,brows,receipt,out,figs,evidence)
    receipt['artifact_generation_seconds']=time.perf_counter()-start-receipt['wall_seconds']
    receipt['total_wall_seconds']=time.perf_counter()-start
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt),flush=True)

def residual_relation(old,out,S):
    """Reuse sealed scalar physical receipts; no additional physical solves."""
    energies={}
    with (out/'per_case_energy.csv').open() as f:
        for r in csv.DictReader(f):
            if r['quantity']!='Maxwell_current':continue
            key=(r['run'],r['obs_id'])
            if key not in energies:energies[key]=dict(meta={k:r[k] for k in META[:11]},states={})
            z=energies[key]['states']
            if r['stage']=='anchor_to_raw_unclipped':
                z['anchor']=float(r['S_field__before']);z['raw_unclipped']=float(r['S_field__after'])
            if r['stage']=='anchor_to_final_top1':z['final_top1']=float(r['S_field__after'])
    records={r['obs_id']:r for r in json.loads((old/'evaluation/CONFIRMATION_SEALED_MANIFEST.json').read_text())['records']}
    obs_cache={};rows=[]
    for (run,obs),v in energies.items():
        if obs not in obs_cache:
            r=records[obs]
            with np.load(old/'data'/r['public_path'],allow_pickle=False) as z:Y=np.array(z['Y'])
            with np.load(old/'data'/r['truth_path'],allow_pickle=False) as z:jt=np.array(z['j'])
            eta=Y-S@jt;obs_cache[obs]=(squared_norm(Y),squared_norm(eta))
        y2,n2=obs_cache[obs]
        receipt=json.loads((old/'evaluation'/run/f'{obs}.json').read_text())
        cand=receipt.get('candidates',[])
        anchor=next((c for c in cand if c['label']=='anchor'),None)
        raw=next((c for c in cand if c['label']=='raw_1'),None)
        with np.load(old/'evaluation'/run/f'{obs}.npz',allow_pickle=False) as z:
            raw_difference=float(np.max(abs(z['m_proposal_unclipped']-z['m_proposal_clipped'])))
            raw_matches=raw_difference==0.
        for stage,t2 in v['states'].items():
            source=anchor if stage=='anchor' else receipt.get('top1') if stage=='final_top1' else raw if raw_matches else None
            base=dict(v['meta'],stage=stage,true_S_field_error_energy=t2,observation_energy=y2,noise_energy=n2,
                eta_definition='eta = Y - S @ historical sealed GT current',
                identity='Fpred-FGT = (Fpred-Y) + eta; true_error = observed_misfit + noise + signed_cross',
                residual_source='frozen evaluator receipt; no new solve')
            if source is None:
                rows.append(dict(base,status='NOT_RUN_EXACT_UNCLIPPED_RAW_RECEIPT_FLOAT32_BOX_ROUNDOFF_ONLY'if raw_difference<1e-7 else'NOT_RUN_RAW_UNCLIPPED_DIFFERS_FROM_VERIFIED_CLIPPED_CANDIDATE',observed_data_relative_residual='',observed_misfit_energy='',signed_residual_noise_cross='',raw_vs_verified_clipped_Linf=raw_difference,verified_clipped_raw_relative_residual=None if raw is None else raw['relative_residual']))
            else:
                rel=float(source['relative_residual']);mis=rel*rel*y2
                rows.append(dict(base,status='RECORDED',observed_data_relative_residual=rel,observed_misfit_energy=mis,signed_residual_noise_cross=t2-mis-n2,raw_vs_verified_clipped_Linf=raw_difference,verified_clipped_raw_relative_residual=None if raw is None else raw['relative_residual']))
    writecsv(out/'data_residual_relation.csv',rows)
    return dict(rows=len(rows),recorded=sum(r['status']=='RECORDED' for r in rows),NOT_RUN=sum(r['status']!='RECORDED' for r in rows),additional_LU=0)

def plots_and_findings(rows,boots,receipt,out,figs,evidence):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    arms=['v_image','v_physics','j_full','j_null','j_null_current']
    # Absolute pooled means, rather than share-of-error bars.
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,(quantity,projection) in zip(axs.flat,[('material_complex_Green_embedding','total'),('material_complex_Green_embedding','D509'),('Maxwell_current','total'),('Maxwell_current','D509')]):
        for noise,marker in [(0.,'o'),(.001,'s')]:
            for stage,style in [('anchor_to_raw_unclipped','--'),('anchor_to_final_top1','-')]:
                vals=[]
                for arm in arms:
                    rs=[r for r in rows if r['arm']==arm and r['seed']==7 and r['passes']==1 and r['noise']==noise and r['level']=='low' and r['quantity']==quantity and r['projection']==projection and r['stage']==stage]
                    vals.append(sum(r['after_sum'] for r in rs)/max(sum(r['cases'] for r in rs),1))
                ax.plot(arms,vals,marker=marker,linestyle=style,label=f'{"clean" if noise==0 else "noisy"}, {"raw" if "raw_unclipped" in stage else "final"}')
        ax.set_yscale('log');ax.set_title(quantity.replace('_complex_Green_embedding','')+' / '+projection);ax.set_ylabel('absolute squared error mean');ax.tick_params(axis='x',rotation=25);ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=8);fig.savefig(figs/'absolute_raw_final.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    for ax,qty in zip(axs,['material_complex_Green_embedding','Maxwell_current']):
        for stage,color,offset in [('anchor_to_raw_unclipped','C0',-.14),('anchor_to_final_top1','C1',.14)]:
            rr=[next((r for r in boots if r['arm']==arm and r['seed']==7 and r['passes']==1 and r['noise']==.001 and r['level']=='low' and r['quantity']==qty and r['projection']=='D509' and r['stage']==stage),None) for arm in arms]
            y=np.array([r['delta_mean'] if r else np.nan for r in rr]);lo=np.array([r['delta_ci_low'] if r else np.nan for r in rr]);hi=np.array([r['delta_ci_high'] if r else np.nan for r in rr])
            ax.errorbar(np.arange(len(arms))+offset,y,yerr=[y-lo,hi-y],fmt='o',color=color,label='raw' if 'raw_unclipped' in stage else 'final')
        ax.axhline(0,color='k',linewidth=.7);ax.set_xticks(range(len(arms)),arms,rotation=25);ax.set_title(qty.replace('_complex_Green_embedding',''));ax.set_ylabel('D509 absolute energy change vs GN');ax.grid(alpha=.2);ax.legend()
    fig.savefig(figs/'noisy_D509_paired_delta.png',dpi=160);plt.close(fig)
    # Signed update/cross decomposition shows where apparent correction arises.
    fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    for ax,qty in zip(axs,['material_complex_Green_embedding','Maxwell_current']):
        for metric,offset in [('update_mean',-.24),('signed_cross_mean',0),('delta_mean',.24)]:
            vals=[]
            for arm in arms:
                rs=[r for r in rows if r['arm']==arm and r['seed']==7 and r['passes']==1 and r['noise']==.001 and r['level']=='low' and r['quantity']==qty and r['projection']=='D509' and r['stage']=='anchor_to_final_top1']
                vals.append(sum(r[metric]*r['cases'] for r in rs)/max(sum(r['cases'] for r in rs),1))
            ax.bar(np.arange(len(arms))+offset,vals,width=.23,label=metric.replace('_mean',''))
        ax.axhline(0,color='k',linewidth=.7);ax.set_xticks(range(len(arms)),arms,rotation=25);ax.set_title(qty.replace('_complex_Green_embedding',''));ax.set_ylabel('D509 absolute energy');ax.legend(fontsize=8)
    fig.savefig(figs/'signed_cross_term.png',dpi=160);plt.close(fig)
    selected=[r for r in boots if r['seed']==7 and r['passes']==1 and r['level']=='low' and r['projection']=='D509' and r['stage'] in ('anchor_to_raw_unclipped','anchor_to_final_top1') and r['quantity'] in ('material_complex_Green_embedding','Maxwell_current')]
    lines=['# E0 absolute energy audit','',f"Status: **{receipt['status']}**. {receipt['counts'].get('frozen_model_observations',0)}/{receipt['expected_model_observations']} historical model-observations audited; no new training or new R4 test truth.",'',
        'Raw denotes the original unclipped NN material proposal. Actual nonlinear Maxwell currents were recomputed for that material and the frozen final top1. Missing physically feasible top1 entries stay NOT_RUN. Auxiliary-head discrepancies are separately labelled and are not Maxwell-solution current errors.','',
        'Green material projections embed real material vectors in complex right-domain Green bases. They are distinct from real material-Jacobian strong directions. D509 is the orthogonal complement of B67 = span(internal strong19, receiver visible48). Sequential I→O, O→I, numerical intersections, and the intersection remainder are separate CSV columns.','',
        'All energy fields are unnormalized squared norms. Each transition records before, after, delta, update, and signed cross; delta = update + signed cross. Main conclusions must use the absolute errors and paired deltas, not energy shares.','',
        '## Low-amplitude D509 paired changes (seed 7, one pass)','',
        '|arm|noise|quantity|stage|mean ΔE|95% mother bootstrap interval|improved mothers|',
        '|---|---:|---|---|---:|---|---:|']
    for r in selected:
        lines.append(f"|{r['arm']}|{r['noise']}|{r['quantity']}|{r['stage']}|{r['delta_mean']:.6g}|[{r['delta_ci_low']:.6g}, {r['delta_ci_high']:.6g}]|{r['mother_improvement_fraction']:.1%}|")
    lines+=['','Bootstrap is paired within mother, with seeds kept separate; low and stress are separate. Every noise/amplitude/arm/seed group is in grouped_absolute_energy.csv. All 17 planned models, including the two seed-7 cascades, use the same fixed geometry. The real-J diagnostic uses only the declared first two parents per family × ID/OOD, both noise observations.','',
        f"Cost: {receipt['wall_seconds']:.1f} s; {receipt['physical_costcounts'].get('LU',0)} LU, {receipt['physical_costcounts'].get('actual_complex_RHS',0)} complex RHS; one CPU worker, BLAS=1. Max full VIE equation defect {receipt['max_full_VIE_equation_defect']:.3g}.",'',
        'True S-field error and observed residual differ: Fpred−FGT = (Fpred−Y)+η. data_residual_relation.csv records observed misfit, noise energy, and signed residual–noise cross from frozen physical receipts, without extra solves. Clean GN anchors can fit much more accurately than the .001 correction stopping gate; increased S-field error can partly reflect this stopping-accuracy difference. This audit does not identify its cause.','',
        'Evidence: per_case_energy.csv; grouped_absolute_energy.csv; mother_paired_bootstrap.csv; data_residual_relation.csv; real_J_diagnostic.csv; J_DIAGNOSTIC_SELECTION.json; GEOMETRY.json; RECEIPT.json; NOT_RUN.csv if any. Scientific interpretation remains with the parent review.']
    (evidence/'E0_FINDINGS.md').write_text('\n'.join(lines)+'\n')
    (out/'E0_FINDINGS.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
