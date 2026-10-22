"""Final evidence figures from scored summaries; never opens truth or LU arrays.

Usage (main schedules under stage D analysis accounting after both evaluations):
  python tools/plot_p3r2_results.py --output-directory /NEW/figure-directory
Missing inputs yield explicit NOT_RUN panels/manifest entries, never zero values.
No campaign imports, fitting, physics, checkpoint loading, or network activity.
"""
from pathlib import Path
import argparse
import json
import time
import zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT/'research/p3l_r2_risk_posterior_exploration'
LIMIT = 32*1024*1024


class Missing(Exception):
    pass


class Figures:
    def __init__(self, root, destination):
        self.root, self.destination = root, destination
        self.sources, self.manifest, self.current = set(), [], None

    def path(self, relative):
        path = self.root/relative
        # All data access is confined to explicit scored-result/source folders.
        if path.is_symlink() or self.root not in path.resolve().parents:
            raise ValueError('External/symlink data refused: '+str(path))
        if any(p in {'inputs', 'offline_only', 'contexts', 'endpoint_geometry'} for p in path.parts):
            raise ValueError('Truth/context/geometry data access forbidden')
        if not path.is_file():
            raise Missing('Missing '+str(relative))
        relative_name=str(path.relative_to(self.root))
        self.sources.add(relative_name)
        if self.current is not None:self.current.setdefault('sources',set()).add(relative_name)
        return path

    def read(self, relative):
        path = self.path(relative)
        if path.stat().st_size > LIMIT:
            raise ValueError('Summary JSON exceeds bounded read limit: '+str(relative))
        return json.loads(path.read_text())

    def arrays(self, relative, *keys):
        path = self.path(relative)
        with zipfile.ZipFile(path) as zipped:
            for key in keys:
                if key+'.npy' not in zipped.namelist():
                    raise Missing('Missing array '+key+' in '+str(relative))
                if zipped.getinfo(key+'.npy').file_size > LIMIT:
                    raise ValueError('Selected array exceeds bounded read limit: '+key)
        with np.load(path, allow_pickle=False) as data:
            return tuple(data[key].copy() for key in keys)

    def missing(self, ax, name, error):
        ax.text(.5, .5, name+'\nNOT_RUN\n'+str(error), ha='center', va='center',
                transform=ax.transAxes, wrap=True, fontsize=9)
        ax.set_axis_off()
        self.current['missing'].append(dict(panel=name, status='NOT_RUN', reason=str(error)))

    def panel(self, ax, name, callback):
        try:
            callback(ax)
        except (Missing, KeyError, IndexError) as error:
            self.missing(ax, name, error)

    def run(self, name, title, shape, callback, scope):
        self.current = dict(id=name, scope=scope, missing=[], outputs=[])
        fig, axes = plt.subplots(*shape, figsize=(6*shape[1], 4.5*shape[0]), squeeze=False,
                                 constrained_layout=True)
        fig.suptitle(title, fontsize=12)
        try:
            callback(axes)
            for extension in ('png', 'pdf'):
                path = self.destination/(name+'.'+extension)
                fig.savefig(path, dpi=180)
                self.current['outputs'].append(path.name)
            self.current['status'] = 'PARTIAL_NOT_RUN' if self.current['missing'] else 'PLOTTED'
            if len(self.current['missing']) == axes.size:
                self.current['status'] = 'NOT_RUN'
        finally:
            plt.close(fig)
            self.current['sources'] = sorted(self.current.get('sources',set()))
            self.manifest.append(self.current)


def points(ax, x, y, label=None, **kwargs):
    pairs = [(a,b) for a,b in zip(x,y) if a is not None and b is not None and np.isfinite(a) and np.isfinite(b)]
    if not pairs:
        raise Missing('No finite observations for '+str(label))
    xx, yy = zip(*pairs)
    ax.plot(xx, yy, 'o-', label=label, **kwargs)


def validation(f, axes):
    groups = [('direct', ['training', 'training32768']), ('residual', ['residual_training']),
              ('decoder controls', ['residual_controls']), ('local controls', ['residual_local_controls'])]
    for ax, (label, folders) in zip(axes.flat, groups):
        def draw(ax, label=label, folders=folders):
            runs=[]
            for folder in folders:
                try:runs.extend(f.read(f'A/{folder}/RESULTS.json')['results'])
                except Missing as e:f.current['missing'].append(dict(panel=label, status='NOT_RUN', reason=str(e)))
            if not runs:raise Missing('No validation summaries')
            family=sorted({r['kind'] for r in runs})
            colors=plt.get_cmap('tab20')
            for k, kind in enumerate(family):
                seeds=sorted({r['seed'] for r in runs if r['kind']==kind})
                for seed in seeds:
                    rr=sorted((r for r in runs if r['kind']==kind and r['seed']==seed),key=lambda r:r['train_size'])
                    points(ax,[r['train_size'] for r in rr],
                           [r['validation']['complex_roi_mse']*256/576 for r in rr],
                           f'{kind} s{seed}',color=colors(k%20),alpha=.7)
                expected={11,22,33}
                if kind not in ('mean','ridge','anchor','residual_ridge') and not expected <= set(seeds):
                    f.current['missing'].append(dict(panel=label+'/'+kind,status='NOT_RUN',reason='Missing declared seeds '+str(sorted(expected-set(seeds)))))
            sizes=sorted({r['validation_size'] for r in runs})
            ax.set(title=f'{label}; validation parents {sizes}',xlabel='Independent training parents',
                   ylabel='Full-grid complex MSE (ROI MSE × 256/576)')
            ax.set_xscale('log',base=2);ax.legend(fontsize=6,ncol=2)
        f.panel(ax,label,draw)


def risk_forest(f, axes):
    result=f.read('L/evaluation/RESULTS.json')
    means=result['summaries']['fullgrid_MSE']
    names=list(means)
    axes[0,0].figure.set_size_inches(15,max(6,.3*len(names)+1.5))
    positions=np.arange(len(names))
    labels=[n+(' [secondary]' if n.startswith('ablation__') or n=='classical_postprocess' else '') for n in names]
    ax=axes[0,0]
    ax.barh(positions,[means[n]['mean'] for n in names])
    ax.set(yticks=positions,yticklabels=labels,xlabel='Parent mean full-grid complex MSE',title='All256; all recorded methods retained')
    ax.tick_params(axis='y',labelsize=7);ax.invert_yaxis()
    comp=result['secondary_comparisons_vs_probit']['fullgrid_MSE']
    ax=axes[0,1]
    used=[n for n in names if n in comp]
    for i,n in enumerate(used):
        row=comp[n];lo,hi=np.asarray(row['relative_reduction_CI_95'])*100
        value=row['relative_risk_reduction']*100
        ax.plot([lo,hi],[i,i],color='C0');ax.plot(value,i,'o',color='C0')
        # Adjusted interval for the prespecified fused co-primary is explicit.
        if n=='fused' and 'relative_reduction_CI_97_5' in row:
            lo2,hi2=np.asarray(row['relative_reduction_CI_97_5'])*100
            ax.plot([lo2,hi2],[i,i],color='C1',alpha=.5,linewidth=4)
    ax.axvline(0,color='grey',linewidth=.8)
    ax.set(yticks=range(len(used)),yticklabels=[labels[names.index(n)] for n in used],
           xlabel='Relative risk reduction vs probit (%)',title='Paired parent bootstrap 95% CI; fused orange 97.5%')
    ax.tick_params(axis='y',labelsize=7);ax.invert_yaxis()
    f.current['fallback_counts']=result.get('per_method_fallback_counts',{})
    f.current['baseline_qualification']=result.get('primary_baseline_qualification')
    f.current['secondary_seed_methods']=[n for n in names if n.startswith('ablation__')]


def density(f, axes):
    result=f.read('L/evaluation/RESULTS.json')['density']
    for i,name in enumerate(('density','weak_density')):
        row=result[name];ax=axes[0,i]
        ax.bar(['learned','unit prior'],[row['mean_nll'],row['unit_prior_mean_nll']])
        ax.set(title=name+f" (prior fallbacks {row['fallback_count']})",ylabel='Proper coordinate NLL; lower is better')
        ax=axes[1,i];coverage=np.asarray(row['component_coverage90'],float)
        ax.bar(np.arange(len(coverage)),coverage);ax.axhline(.9,color='black',linestyle='--')
        ax.set(xlabel='Frozen coordinate (z8 + probit u2)',ylabel='Empirical marginal 90% coverage',ylim=(0,1.02))


def source_smc(f, axes):
    for parent in ('000','030'):
        rel=f'B/jobs/B_source_DEV_{parent}/LANDSCAPE_PROGRESS.json'
        try:
            rr=f.read(rel)['rows']
            if not rr:raise Missing('Empty source SMC progress')
            x=[r['stage'] for r in rr]
            points(axes[0,0],x,[r['distinct_prior_ancestors'] for r in rr],'DEV_'+parent)
            points(axes[0,1],x,[sum(r['target_coefficients']) for r in rr],'DEV_'+parent)
            points(axes[1,0],x,[r.get('minimum_full_energy_excess_over_paid_anchor') for r in rr],'DEV_'+parent)
            points(axes[1,1],[r['new_RHS'] for r in rr],[r['distinct_prior_ancestors'] for r in rr],'DEV_'+parent)
        except (Missing,KeyError) as e:
            f.current['missing'].append(dict(panel='source '+parent,status='NOT_RUN',reason=str(e)))
    for ax,title,xlabel,ylabel in zip(axes.flat,
            ('Prior ancestry','Source/beta progression','Energy excess over paid anchor','Ancestry vs paid work'),
            ('stage','stage','stage','new RHS'),
            ('Distinct prior ancestors','Sum of eight source powers (FULL = 8)','Minimum FULL Phi − paid anchor Phi','Distinct prior ancestors')):
        ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
        if ax.lines:ax.legend()
        else:f.missing(ax,title,'No source progress rows')
    f.current['scope_note']='Sum of source powers records sequential beta progression; not one global inverse temperature. Full energy is descriptive even at partial target.'


def mala(f, axes):
    result=f.read('B/jobs/B_metric_pair_DEV030/RESULT.json')
    arms=result['arms']
    for name,row in arms.items():
        chains=row['chains']
        points(axes[0,0],range(len(chains)),[c['production_acceptance'] for c in chains],name)
        checkpoints=result.get('matched_production_RHS_checkpoints',[])
        for ax,key in ((axes[0,1],'rhat_maximum_finite'),(axes[1,0],'complex_material_mean_RMS_MCSE')):
            x=[];y=[]
            for checkpoint in checkpoints:
                if name in checkpoint['arms']:
                    item=checkpoint['arms'][name];x.append(item['actual_production_RHS_per_chain']);y.append(item['material_functions'].get(key))
            if x:
                try:points(ax,x,y,name)
                except Missing as error:
                    f.current['missing'].append(dict(panel=name+'/'+key,status='NOT_RUN',reason=str(error)))
            else:f.current['missing'].append(dict(panel=name+'/'+key,status='NOT_RUN',reason='No matched-production checkpoints'))
        axes[1,1].bar(name,row['actual_RHS'])
    for ax,title,xlabel,ylabel in zip(axes.flat,
        ('Acceptance','Material functions: maximum finite Rhat','Asymptotic material mean MCSE','Total actual arm cost'),
        ('chain','production RHS per chain','production RHS per chain','arm'),
        ('Production acceptance','Max finite Rhat','Complex RMS MCSE','Actual RHS including warmup')):
        ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
        if ax.lines:ax.legend()
    f.current['undefined_and_stuck']={n:r['material_functions'] for n,r in arms.items()}


def slices(f, axes):
    jobs=['B_field2_slice_DEV030_r0','B_partition_DEV_000','B_partition_DEV_030']
    for name in jobs:
        def draw(ax,name=name):
            result=f.read(f'B/jobs/{name}/RESULT.json');rr=result['quadrature_rows']
            if not rr:raise Missing('No quadrature rows')
            logkey='logZ' if 'logZ' in rr[0] else 'log_integral'
            points(ax,[r['order'] for r in rr],[r[logkey] for r in rr],'log integral')
            ax.set(title=name+'\nUNQUALIFIED conditional SLICE',xlabel='Quadrature order',ylabel='log integral')
            twin=ax.twinx()
            try:
                points(twin,[r['order'] for r in rr],[r.get('material_mean_successive_RMS') for r in rr],
                       'Successive mean RMS',color='C1')
            except Missing as error:
                twin.text(.5,.5,'Successive RMS: NOT_RUN',transform=twin.transAxes,ha='center',fontsize=8)
                f.current['missing'].append(dict(panel=name+'/successive_RMS',status='NOT_RUN',reason=str(error)))
            twin.set_ylabel('Successive material mean RMS',color='C1')
            f.current.setdefault('qualification',{})[name]=dict(reference=result.get('reference_qualification'),stopping_rule_met=result.get('quadrature_stopping_rule_met'))
        f.panel(axes.flat[jobs.index(name)],name,draw)
    def shape(ax):
        job='B/jobs/B_field2_slice_DEV030_r0';result=f.read(job+'/RESULT.json');order=result['final_GH_order']
        a,phi=f.arrays(f'{job}/GH_{order}_RAW.npz','a','Phi')
        value=-phi-.5*np.sum(a*a,axis=1);value-=np.max(value)
        image=ax.scatter(a[:,0],a[:,1],c=value,s=12,cmap='viridis')
        ax.figure.colorbar(image,ax=ax,label='Log density relative to sampled maximum')
        ax.set(title='Field2 conditional density at GH nodes; UNQUALIFIED',xlabel='Gaussian coordinate a1',ylabel='Gaussian coordinate a2')
    f.panel(axes[1,1],'field2 density shape',shape)


def jensen(f, axes):
    job='B/jobs/B_field2_slice_DEV030_r0';result=f.read(job+'/RESULT.json');order=result['final_GH_order']
    gap,variance=f.arrays(f'{job}/GH_{order}_MOMENTS.npz','decoder_Jensen_gap','variance')
    if gap.shape!=(576,) or variance.shape!=(576,):raise ValueError('Expected 24×24 field2 moment maps')
    for ax,data,title in zip(axes.flat,(gap.real,gap.imag,variance),('Real E[D(q)] − D(E[q])','Imag E[D(q)] − D(E[q])','Conditional complex material variance')):
        if 'variance' in title:image=ax.imshow(data.reshape(24,24),origin='lower',cmap='viridis',vmin=0)
        else:
            limit=float(np.max(np.abs(data)))
            image=ax.imshow(data.reshape(24,24),origin='lower',cmap='RdBu_r',vmin=-limit if limit else -1e-15,vmax=limit if limit else 1e-15)
        ax.set_title(title);ax.set_xticks([]);ax.set_yticks([]);ax.figure.colorbar(image,ax=ax,shrink=.8)
    f.current['reference_qualification']=result.get('reference_qualification')


def trajectories(f, axes):
    result=f.read('L/trajectory_evaluation/RESULTS.json')
    core=result['primary_arm'];risk=result['per_parent_risks'];parents=result['parents']
    for ax,cohort,label in ((axes[0,0],parents,'all256 core arm'),(axes[0,1],parents[:32],'same first32 secondary arms/controls')):
        arms=[core] if len(cohort)==256 else list(result['summaries'])
        base=np.mean([risk[p]['probit']['fullgrid_MSE'] for p in cohort])
        ax.axhline(base,label='probit risk (cost not assigned)',color='grey',linestyle='--')
        for arm in arms:
            diag=result['actual_candidate_diagnostics'][arm]
            # Cost is recorded for the full registered arm cohort, not assumed
            # identical across its parents or across core/control arms.
            if len(cohort)==256 or arm!=core:
                x=diag['incremental_step_RHS']/diag['parents']
            else:
                # All256 cost cannot be attributed to the first32 without
                # per-parent endpoint metadata. No fabricated proportional cost.
                statuses=[]
                for p in cohort:
                    _,=f.arrays(f'L/trajectory_estimates/{p}.npz','status_json')
                    statuses.append(json.loads(str(_)))
                x=sum(sum(r['cost']['linear_solve_rhs'] for r in s['arms'][arm]['selected_endpoint_diagnostics']) for s in statuses)/len(cohort)
            for i,method in enumerate(result['summaries'][arm]):
                name=arm+'__'+method;y=np.mean([risk[p][name]['fullgrid_MSE'] for p in cohort])
                ax.scatter(x,y,marker=('o','s','^','D','v','P')[i%6],label=arm+'/'+method)
        ax.set(title=label,xlabel='Recorded incremental trajectory RHS / parent\n(context, setup and training costs excluded)',ylabel='Full-grid complex MSE')
        ax.legend(fontsize=6)
    f.current['cost_scope']='Incremental selected-step RHS only. No equal-total-cost or posterior integration claim.'


def forwards(f, axes):
    result=f.read('L/jobs/locked_decision_forward/RESULTS.json')
    if result.get('status')!='COMPLETE':raise Missing('Decision-forward job not COMPLETE')
    rr=result['rows'];methods=sorted({r['method'] for r in rr})
    for name in methods:
        values=np.asarray([r['Phi'] for r in rr if r['method']==name and r.get('status')=='COMPLETE' and 'Phi' in r],float)
        if not values.size:
            f.current['missing'].append(dict(panel=name,status='NOT_RUN',reason='No completed forward values'));continue
        values=np.sort(values)
        axes[0,0].step(values,np.arange(1,len(values)+1)/len(values),where='post',label=f'{name} n={len(values)}')
        axes[0,1].hist(values,bins=min(30,max(1,int(np.sqrt(len(values))))),histtype='step',label=name)
    for ax in axes.flat:
        ax.set_xlabel('Phi of saved deterministic risk action');ax.legend(fontsize=7)
    axes[0,0].set_ylabel('Empirical CDF');axes[0,1].set_ylabel('Count')
    f.current['unevaluated_rows']=[r for r in rr if r.get('status')!='COMPLETE']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--output-directory',required=True,type=Path,help='NEW directory, never overwrite final figures')
    args=parser.parse_args();root=args.root.resolve();dest=args.output_directory.expanduser().resolve()
    if dest.exists():raise FileExistsError('Choose a NEW figure directory')
    f=Figures(root,dest)
    markers=[f.read('L/'+name+'/COMPLETE.json') for name in ('evaluation','trajectory_evaluation')]
    for m in markers:
        if m.get('status')!='COMPLETE' or m.get('count')!=256 or len(m.get('parents',[]))!=256 or len(set(m.get('parents',[])))!=256:
            raise ValueError('Both exact256 evaluation COMPLETE gates required')
    for key in ('parents','source_revision','frozen_unix'):
        if key not in markers[0] or markers[0][key]!=markers[1].get(key):raise ValueError('Evaluation freeze mismatch: '+key)
    dest.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    specifications=[
        ('01_validation','Validation learning curves; individual seeds, distinct validation cohorts',(2,2),validation,'ROI-to-fullgrid factor 256/576; no full-posterior inference'),
        ('02_locked_risk','Fresh locked risk; all256 parents retained',(1,2),risk_forest,'95% paired intervals; all ablations/seed controls secondary; no test reselection'),
        ('03_coordinate_density','Frozen-coordinate proper scores and marginal coverage only',(2,2),density,'Not FULL calibration; unit Gaussian coordinate prior baseline'),
        ('04_source_smc','Source SMC ancestry and energy; explicit partial targets retained',(2,2),source_smc,'No posterior-reference qualification; energy excess is descriptive'),
        ('05_metric_mala','Frozen/dynamic Full-J MALA diagnostics; UNQUALIFIED',(2,2),mala,'Finite Rhat/MCSE summaries retain undefined functions; acceptance is not mixing'),
        ('06_slice_convergence','Numerical conditional SLICE refinement; UNQUALIFIED',(2,2),slices,'Field2 Gaussian and property U2 conditional slices; not FULL target reference'),
        ('07_slice_jensen','Field2 conditional SLICE moments; UNQUALIFIED, not FULL',(1,3),jensen,'Decoder Jensen gap and conditional variance at final GH order'),
        ('08_trajectory_risk','Fresh trajectory risk versus recorded incremental work',(1,2),trajectories,'All256 core; identical first32 secondary panel; coefficients are estimator weights'),
        ('09_decision_forward','Forward of deterministic risk actions; not PPC',(1,2),forwards,'No posterior predictive or decoder legality claim')]
    complete=False
    try:
        for name,title,shape,callback,scope in specifications:
            def guarded(axes,callback=callback):
                try:callback(f,axes)
                except (Missing,KeyError,IndexError) as error:
                    for ax in axes.flat:f.missing(ax,name,error)
            f.run(name,title,shape,guarded,scope)
        complete=True
    finally:
        manifest=dict(status='FIGURES_WRITTEN_WITH_EXPLICIT_MISSING_DATA' if complete else 'FAILED_PARTIAL_OUTPUT_RETAINED',truth_file_access=False,
                      physics_calls=0,source_revision=markers[0]['source_revision'],frozen_unix=markers[0]['frozen_unix'],
                      wall_seconds=time.monotonic()-start,figures=f.manifest,all_sources=sorted(f.sources),
                      limitations='Source implementation is not scientific acceptance; PNG/PDF rendering requires visual review.')
        with (dest/'FIGURE_MANIFEST.json').open('x') as stream:json.dump(manifest,stream,indent=2,allow_nan=False)
    print(json.dumps(dict(directory=str(dest),figures=len(f.manifest),manifest='FIGURE_MANIFEST.json')))


if __name__=='__main__':main()
