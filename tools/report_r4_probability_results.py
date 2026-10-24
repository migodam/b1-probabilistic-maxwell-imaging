"""Read CLOSED scores only: probability tables/figures, no new score or physics."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import csv,json,resource,shutil,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from scipy.special import gammainc
from scipy.stats import gamma
from p3r4.core import Book,OUT,write
from p3r4.train_physics_transport import identity,verify_identity

GATE=OUT/'results/FRESH_ISOLATED_SCORES_COMPLETE.json'
TABLES=OUT/'AUDIT/final_probability'
FIGURES=OUT/'figures/final_probability'
LABELS={'raw_spline':'Raw spline','matched_diag':'Matched diagonal Gaussian','tri_maxwell':'Triangular TRAIN Fisher frame',
 'tri_random':'Triangular fixed random frame','raw_cfm__nfe32':'Raw CFM','anchor_cfm_none__nfe32':'Anchor CFM',
 'anchor_cfm_geometry__nfe32':'Geometry-weighted CFM','anchor_cfm_random__nfe32':'Random-alignment CFM',
 'refined_physics':'Physics-refined spline (S1)','refined_supervised':'Supervised-continuation spline (S1)'}
COLORS=['#2166ac','#d95f02','#1b9e77','#7b3294']


def load(path):return json.loads(Path(path).read_text())


def csvout(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (TABLES/name).open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=keys);w.writeheader();w.writerows(rows)


def savefig(fig,name):
    fig.savefig(FIGURES/(name+'.png'),dpi=200,bbox_inches='tight')
    fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight');plt.close(fig)


def mean_metric(container,key):return container.get('metrics',{}).get(key,{}).get('mean_all_parents')


def errorbar(summary):
    mean=summary.get('mean_all_parents');ci=summary.get('CI')
    return np.array([[mean-ci[0]],[ci[1]-mean]]) if mean is not None and ci else None


def closed_gate():
    if load(GATE).get('status')!='COMPLETE':raise ValueError('global isolated scoring barrier not COMPLETE')
    release=OUT/'configs/FIRST_HELDOUT_GT_SCORE_RELEASE.json';r=load(release)
    if r.get('released') is not True or r.get('experiment')!='fresh_r4_isolated_scores_v1':raise ValueError('first heldout GT release missing')
    for record in r['closed_gates']:verify_identity(record)
    receipts={}
    for name in ('fresh_density_locked_score_NFE32','fresh_density_confirmation_score_NFE32','fresh_full_probability_offline_score'):
        p=OUT/'IV/jobs'/name/'CLOSED.json'
        if load(p)['status']!='COMPLETE':raise ValueError('respective CLOSED score gate required: '+name)
        receipts[name]=identity(p)
    return r,receipts


def metadata(frozen):
    import torch
    top=load(OUT/'configs/FROZEN_DENSITY_BOTH_COHORTS.json')
    originals={m['id']:m for m in top['methods']};rows={}
    for m in frozen['methods']:
        verify_identity(m['checkpoint_identity']);verify_identity(originals[m['id']]['source_receipt'])
        p=torch.load(m['checkpoint'],map_location='cpu',weights_only=False);c=p['contract'];a=p.get('anchor_contract',{})
        family=m['family'];frame=c.get('frame',c.get('frame_path'))
        description='fixed original coordinate order'
        if family=='tri_maxwell':description='TRAIN256 likelihood-Fisher ascending full66 frame; first16 conditional block; no per-test anchor'
        elif family=='tri_random':description='fixed TRAIN-independent orthogonal frame seed20261025113; first16 conditional block'
        elif family=='anchor_cfm_geometry':description='coordinate none; TRAIN standardized Jr=Jx*T bounded SPD velocity-error weighting; no test-J feature'
        elif family=='anchor_cfm_random':description='coordinate none; same per-parent metric spectrum, fixed conjugation R seed20261025191; no test-J feature'
        elif family.startswith('anchor_cfm'):description='coordinate none; raw512+observation-derived anchor66; fullrank computational base'
        rows[m['id']]=dict(id=m['id'],family=family,seed=m['seed'],kind=m['kind'],
            TRAIN=c['train_size'],VAL=c['validation_size'],condition_dim=p['model_config']['context_dim'],
            parameter_count=c['parameter_count'],selected_update=p['selected_update'],selection_criterion=p['criterion'],
            fixed_frame=frame,frame_information=description,
            physics_TRAIN_pool=c.get('physics_train_size',0),
            validation_reverse_roster='32 parents x4 draws' if 'refinement_contract' in p else 'NOT_APPLICABLE',
            reference_measure='raw Gaussianized dx66',checkpoint=m['checkpoint'])
        del p
    return rows


def tables_and_bars(cohorts,meta):
    modelrows=[];familyrows=[];paired=[]
    for cohort,r in cohorts.items():
        for mid,m in r['methods'].items():
            original=mid.removesuffix('__nfe32');row=dict(cohort=cohort,parents=r['n_parent'],variant=mid,**meta[original])
            row.update(NLL_status=m['full66_nll_status'],proper_NLL=mean_metric(m,'full66_nll'),
                material_energy=mean_metric(m,'material_energy_score'),samplemean_MSE=mean_metric(m,'raw_sample_mean_full_mse'),
                phase_Brier=mean_metric(m,'phase_brier_per_threshold'),successful_parents=m['successful_parents'],failures=m['failures'])
            if m['kind']=='cfm':
                if row['proper_NLL'] is not None:raise ValueError('CFM numeric NLL must not be presented as proper exact NLL')
                row['NLL_status']='NOT_RUN_UNQUALIFIED_CFM_NUMERIC_DENSITY'
            for level in (50,80,90,95):
                row[f'ROI_empirical_coverage_{level}']=mean_metric(m,f'material_real_roi__empirical_interval_coverage_{level}')
                row[f'ROI_randomized_rank_coverage_{level}']=mean_metric(m,f'material_real_roi__randomized_rank_coverage_{level}')
            modelrows.append(row)
        for family,data in r['families'].items():
            row=dict(cohort=cohort,parents=r['n_parent'],family=family,seeds=';'.join(map(str,data['training_seeds'])),
                aggregation='per-parent seed-loss mean',proper_NLL=mean_metric(data,'full66_nll'),
                NLL_status='NOT_RUN' if family.endswith('__nfe32') else 'EVALUATED',
                material_energy=mean_metric(data,'material_energy_score'),samplemean_MSE=mean_metric(data,'raw_sample_mean_full_mse'),
                phase_Brier=mean_metric(data,'phase_brier_per_threshold'))
            for key,value in data['metrics'].items():row[key]=value.get('mean_all_parents')
            familyrows.append(row)
        for comp in r['comparisons']:
            item=dict(cohort=cohort,**comp)
            if comp['status']=='EVALUATED' and comp['metric']=='full66_nll':
                item['expected_forward_KL_reduction_nats_estimate']=comp['improvement_absolute']
            if comp['status']=='EVALUATED' and comp['metric']=='material_energy_score':
                item['expected_conditional_energy_distance_reduction_estimate']=2*comp['improvement_absolute']
                item['interpretation']='absolute difference only; score-percentage is not unknown posterior-distance-percentage'
            paired.append(item)
    csvout('MODEL26_BY_COHORT.csv',modelrows);csvout('FAMILY_SEEDLOSS_BY_COHORT.csv',familyrows)
    csvout('FROZEN_PAIRED_COMPARISONS.csv',paired);csvout('FRAME_AND_DATA_INFO.csv',list(meta.values()))
    write(TABLES/'TABLES.json',dict(models=modelrows,families=familyrows,comparisons=paired))
    order=list(next(iter(cohorts.values()))['families'])
    fig,axes=plt.subplots(1,2,figsize=(13,7),sharey=True)
    eligible=[mean_metric(data,'full66_nll') for r in cohorts.values() for data in r['families'].values() if mean_metric(data,'full66_nll') is not None]
    lower=min(0.,min(eligible)-5);upper=max(eligible)*1.15
    for ax,(cohort,r) in zip(axes,cohorts.items()):
        for j,family in enumerate(order):
            summary=r['families'][family]['metrics'].get('full66_nll',{})
            value=summary.get('mean_all_parents')
            if value is None:ax.text(lower+.57*(upper-lower),j,'NOT_RUN (CFM)',va='center',color='#777777',fontsize=8)
            else:ax.barh(j,value,xerr=errorbar(summary),height=.62,color='#4477aa',capsize=2)
        ax.set(xlim=(lower,upper),xlabel='Proper full66 NLL (raw dx66, nats)',title=f'{cohort}: n={r["n_parent"]} parents')
        ax.grid(axis='x',alpha=.2);ax.set_yticks(range(len(order)),[LABELS.get(x,x) for x in order])
    axes[0].invert_yaxis();fig.suptitle('Frozen full66 densities: conditional proper log score; CFM unavailable, not zero')
    fig.tight_layout();savefig(fig,'proper_NLL_status')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    arms=['anchor_cfm_none__nfe32','anchor_cfm_geometry__nfe32','anchor_cfm_random__nfe32']
    for ax,(cohort,r) in zip(axes,cohorts.items()):
        for j,family in enumerate(arms):
            s=r['families'][family]['metrics']['material_energy_score'];value=s['mean_all_parents']
            if value is not None:ax.bar(j,value,yerr=errorbar(s),color=COLORS[j],capsize=3)
        ax.set_xticks(range(3),['Baseline','Geometry','Random']);ax.set(ylabel='Material energy score',title=f'{cohort}, parent n={r["n_parent"]}')
        ax.grid(axis='y',alpha=.2)
    fig.suptitle('K128, NFE32; three seed losses averaged within parent')
    fig.tight_layout();savefig(fig,'material_proper_energy_CFM_controls')
    return modelrows,familyrows,paired


def sample_panels():
    cohorts=OUT/'results/fresh/locked/density_online'
    for parent in ['R4L1_00000','R4L1_00001']:
        with np.load(OUT/'inputs/offline_only'/(parent+'.npz'),allow_pickle=False) as z:truth=z['chi'].copy()
        draws=[]
        for model in ('raw_spline_s2026102501','anchor_cfm_geometry_s2026102501__nfe32'):
            with np.load(cohorts/model/(parent+'.npz'),allow_pickle=False) as z:
                if not np.array_equal(z['sample_index'],np.arange(128)):raise ValueError('fixed original sample order required')
                draws.append(dict(samples=z['chi_samples'][:4],mean=z['material_mean'],vr=z['material_variance_real'],vi=z['material_variance_imag']))
        for component,minimum,maximum,vmax,threshold in [('Re',.1,1.9,.81,1e-4),('Im',.04,.22,.0081,1e-6)]:
            extract=np.real if component=='Re' else np.imag
            fig,axes=plt.subplots(2,7,figsize=(17,5))
            artists=[];variance_artists=[]
            for row,(label,d) in enumerate(zip(['Raw spline S1','Geometry CFM S1'],draws)):
                fields=[truth,*d['samples'],d['mean']]
                for col,array in enumerate(fields):
                    artists.append(axes[row,col].imshow(extract(array).reshape(24,24),origin='lower',cmap='viridis',vmin=minimum,vmax=maximum,interpolation='nearest'))
                variance=d['vr'] if component=='Re' else d['vi']
                variance_artists.append(axes[row,6].imshow(variance.reshape(24,24),origin='lower',cmap='magma',
                    norm=SymLogNorm(linthresh=threshold,vmin=0,vmax=vmax),interpolation='nearest'))
                axes[row,0].set_ylabel(label,fontsize=10)
                for ax in axes[row]:ax.set_xticks([]);ax.set_yticks([])
            for ax,title in zip(axes[0],['GT','Fixed sample 0','Fixed sample 1','Fixed sample 2','Fixed sample 3','K128 mean','Variance (ddof0)']):ax.set_title(title,fontsize=10)
            fig.colorbar(artists[0],ax=axes[:,:6],location='bottom',fraction=.08,pad=.06,shrink=.7,label=f'{component} chi; shared fixed range')
            fig.colorbar(variance_artists[0],ax=axes[:,6],location='bottom',fraction=.08,pad=.06,label=f'Var({component} chi), symlog')
            fig.suptitle(f'{parent}: fixed q-model samples (UNQUALIFIED); no GT-based selection',fontsize=13)
            savefig(fig,parent+'_'+component+'_fixed_samples')


def full_calibration():
    root=OUT/'results/fresh/full_probability';r=load(root/'scores/RESULTS.json');f=load(root/'FROZEN.json')
    groups=f['generated'];phis={};ranks={};rows=[]
    rank_lookup={(a['group'],a['parent']):a['data_dependent_likelihood_rank'] for a in r['records']}
    for group in groups:
        values=[];rr=[]
        for parent in group['parents']:
            with np.load(root/'physics'/group['id']/parent/'LIKELIHOOD_VECTOR.npz',allow_pickle=False) as z:
                indices=z['sample_indices'];phi=z['phi']
                if indices[0]!=0:raise ValueError('fixed0 must remain first')
                values.append(float(phi[0]))
            rank=rank_lookup[(group['id'],parent)]
            if rank.get('status')=='EVALUATED':rr.append(rank['randomized_inclusive_rank'])
        phis[group['id']]=np.asarray(values);ranks[group['id']]=np.asarray(rr)
        rows.append(dict(group=group['id'],**r['Gamma256'][group['id']],likelihood_rank_evaluated=len(rr),
            likelihood_rank_planned=16 if group['id'].startswith('locked__') else 0,
            normalized_importance_density=group['density'],posterior_qualification='NOT_ESTABLISHED'))
    csvout('FULL_CALIBRATION_FIXED0_AND_SBC.csv',rows);write(TABLES/'FULL_CALIBRATION.json',dict(rows=rows,ranks={k:v.tolist() for k,v in ranks.items()}))
    positive=np.concatenate([v[np.isfinite(v)&(v>0)] for v in phis.values()]);lo=min(1.,np.floor(np.log10(positive.min())));hi=max(3.,np.ceil(np.log10(positive.max())))
    bins=np.linspace(lo,hi,45);grid=np.linspace(lo,hi,1000);curve=gamma.pdf(10**grid,a=256)*(10**grid)*np.log(10)
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for row,cohort in enumerate(('locked','confirmation')):
        chosen=[g for g in groups if g['id'].startswith(cohort+'__')]
        for j,g in enumerate(chosen):
            a=phis[g['id']];finite=np.isfinite(a);label=g['id'].split('__',1)[1].replace('_s2026102501','').replace('__nfe32','')
            axes[row,0].hist(np.log10(a[finite&(a>0)]),bins=bins,density=True,histtype='step',linewidth=1.8,color=COLORS[j],label=f'{label}: {finite.sum()}/{len(a)}')
            axes[row,1].hist(gammainc(256,a[finite]),bins=np.linspace(0,1,21),density=True,histtype='step',linewidth=1.8,color=COLORS[j],label=label)
        axes[row,0].plot(grid,curve,'k--',label='Gamma256 reference');axes[row,1].axhline(1,color='k',linestyle='--')
        axes[row,0].set(title=f'{cohort}: fixed sample0 Full Phi',xlabel='log10 Phi',ylabel='Density in log10 Phi')
        axes[row,1].set(title=f'{cohort}: Gamma256 PIT',xlabel='Gamma256 CDF(Phi)',ylabel='Density')
        for ax in axes[row]:ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Frozen necessary predictive checks; failures retained, no Full posterior certificate')
    fig.tight_layout();savefig(fig,'fixed0_FullPhi_Gamma256_and_PIT')
    fig,ax=plt.subplots(figsize=(10,4))
    for j,g in enumerate([g for g in groups if g['id'].startswith('locked__')]):
        a=ranks[g['id']]
        if len(a):ax.hist(a,bins=np.linspace(0,1,11),density=True,histtype='step',linewidth=2,color=COLORS[j],label=f'{g["id"].split("__",1)[1]} ({len(a)}/16)')
    ax.axhline(1,color='k',linestyle='--');ax.set(xlabel='Randomized inclusive likelihood rank',ylabel='Density',
        title='Fixed first16 L1 parents, first64 draws: data-dependent SBC diagnostic')
    ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();savefig(fig,'fixed_likelihood_SBC_L1_16')


def vi_plots():
    registry=load(OUT/'configs/PHYSICS_TRANSPORT_PAIRED_V2_RUNTIME_RECOVERY.json')
    families=registry['families'];fig,axes=plt.subplots(3,3,figsize=(13,10))
    scatter,axs=plt.subplots(1,3,figsize=(13,4))
    for row,family in enumerate(families):
        for mode,color in [('physics',COLORS[3]),('supervised',COLORS[0])]:
            name=f"{registry.get('job_prefix','transport')}_{family['id']}_{mode}_s2026102501"
            receipt=OUT/'I/jobs'/name/'CLOSED.json'
            if load(receipt)['status']!='COMPLETE':raise ValueError('VI curves must be CLOSED')
            r=load(Path(registry['output_root'])/name/'RESULT.json');h=r['history'];updates=[v['update'] for v in h]
            for col,(key,label) in enumerate([('validation_full_NLL','Full proper NLL'),('validation_Phi','Mean Full Phi'),('validation_entropy','MC entropy (raw66 nats)')]):
                values=[v[key] for v in h];axes[row,col].plot(updates,values,'o-',ms=3,color=color,label=mode)
                chosen=next(v for v in h if v['update']==r['selected_update'])
                axes[row,col].plot(chosen['update'],chosen[key],'*',color=color,ms=12)
                axes[row,col].set(xlabel='Update',ylabel=label,title=f'{family["id"]}: TRAIN{family["train_size"]}/VAL{family["validation_size"]}')
                axes[row,col].grid(alpha=.2)
            axs[row].plot([v['validation_full_NLL'] for v in h],[v['validation_Phi'] for v in h],'o-',ms=3,color=color,label=mode)
            axs[row].set(xlabel='Proper VAL NLL (nats)',ylabel='Fixed-roster VAL Full Phi',title=family['id']);axs[row].grid(alpha=.2)
        axes[row,0].legend();axs[row].legend()
    fig.suptitle('Closed paired refinement: frozen composite selection (stars); all three replication gates failed')
    fig.tight_layout();savefig(fig,'paired_VI_NLL_Phi_entropy')
    scatter.suptitle('Validation tradeoffs, same fixed128 draws; panels have distinct proper-NLL VAL sizes')
    scatter.tight_layout();savefig(scatter,'paired_VI_NLL_vs_Phi')


def chronology(release,cohorts,score_receipts):
    first=load(OUT/'results/fresh/scoring_queue_v1/FIRST_HELDOUT_TRUTH_GATE.json');t=first['time_ns']
    entries=[]
    for cohort,r in cohorts.items():
        online=load(OUT/f'results/fresh/{cohort}/density_online/ONLINE_COMPLETE.json')
        f=load(OUT/f'results/fresh/{cohort}/FROZEN_DENSITY.json')
        observations=[(OUT/'inputs/online'/(p+'.npz')).stat().st_mtime_ns for p in f['parents']]
        entries.append(dict(cohort=cohort,n_parent=len(f['parents']),density_freeze_time_ns=f['freeze_time_ns'],
            earliest_observation_mtime_ns=min(observations),all_observations_after_density_freeze=min(observations)>=f['freeze_time_ns'],
            online_complete_time_ns=online['completed_time_ns'],online_before_first_GT=online['completed_time_ns']<t,
            inference_pid=online['pid'],score_pid=r['scoring_pid'],distinct_infer_score_PID=online['pid']!=r['scoring_pid'],
            first_GT_score_gate_time_ns=t))
    result=dict(status='CHECKED_RECORDED_CHRONOLOGY',first_GT_release=identity(OUT/'configs/FIRST_HELDOUT_GT_SCORE_RELEASE.json'),
        global_score_gate=identity(GATE),all_release_closed_gates_before_first_GT=all(s['mtime_ns']<=t for s in release['closed_gates']),
        density=entries,respective_closed_score_receipts=score_receipts,model_order='frozen order; no GT reranking',
        sample_panels='first2 L1 objects, fixed indices0..3, q-model mean/variance; no candidate selection')
    write(TABLES/'CHRONOLOGY.json',result);return result


def main():
    while not GATE.exists():time.sleep(60)
    release,score_receipts=closed_gate()
    note=OUT/'AUDIT/JOINT_SIMULATION_PROPER_SCORE_INTERPRETATION.md'
    interpretation=note.read_text()
    if TABLES.exists() or FIGURES.exists():raise FileExistsError('new final probability artifacts only; no overwrite')
    TABLES.mkdir(parents=True);FIGURES.mkdir(parents=True)
    job='final_probability_results_'+str(time.time_ns());book=Book('IV',job,analysis_only=True,max_wall=900)
    started=time.monotonic();status,error='COMPLETE',None
    try:
        cohorts={c:load(OUT/f'results/fresh/{c}/density_scores/RESULTS.json') for c in ('locked','confirmation')}
        if any(r['status']!='OFFLINE_SCORE_COMPLETE' for r in cohorts.values()):raise ValueError('closed isolated density scores required')
        frozen=load(OUT/'results/fresh/locked/FROZEN_DENSITY.json');meta=metadata(frozen)
        modelrows,familyrows,paired=tables_and_bars(cohorts,meta)
        sample_panels();full_calibration();vi_plots();chron=chronology(release,cohorts,score_receipts)
        shutil.copyfile(note,TABLES/'PARENT_PROPER_SCORE_INTERPRETATION.md')
        adaptive=OUT/'AUDIT/adaptive_cfm_density_1792881444063329000'
        for filename in ('REPORT.md','RESULT.json','CONVERGENCE.csv','SOURCE_VERIFICATION.json'):
            shutil.copyfile(adaptive/filename,TABLES/('ADAPTIVE_CFM_'+filename))
        top=load(OUT/'configs/FROZEN_DENSITY_BOTH_COHORTS.json')
        for m in top['methods']:verify_identity(m['checkpoint_identity'])
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>2*1024**3:raise MemoryError('2GiB reporting cap')
        write(TABLES/'RESULT.json',dict(status='COMPLETE',models_per_cohort=26,cohorts=[512,256],model_rows=len(modelrows),
            family_rows=len(familyrows),score_comparison_interpretation='relative expected forward-KL or conditional material energy-distance differences; no absolute posterior-distance percentages',
            primary_source_interpretation=identity(note),chronology=chron,models_unchanged=True,GT_reranking=False,
            generated_new_samples=0,physics_RHS=0,training_updates=0,duplicate_scoring=False,posterior_qualification='NOT_ESTABLISHED',
            wall_seconds=time.monotonic()-started,receipt=str(book.jobdir/'CLOSED.json'),figures=str(FIGURES)))
        print(TABLES,FIGURES,flush=True)
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
