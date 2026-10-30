#!/usr/bin/env python3
"""Separate post-run DEVELOPMENT scorer; online predictions never receive truth."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]))
from research.research_z_alg_r3_digits.experiments.score_predictions import metrics
HISTORY=ROOT.parent/'research_z_alg_r3_digits'
METHODS=['trf','data_lbfgs','globaldual','basin_bridge','ray_best_only']
COLORS=dict(zip(METHODS,['#0072b2','#009e73','#d55e00','#cc79a7','#e69f00']))
COLORS.update(uniform_bridge='#882255',uniform_best_only='#44aabb',uniform_multistart='#332288')

def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def score_run(path,factor):
    manifest=json.loads((path/'MANIFEST.json').read_text())
    rows=[];materials={};runs={};truths={}
    mapping={r['case_id']:r for r in json.loads((ROOT/'truth_development/mapping.json').read_text())['records']} if factor else {}
    for job in manifest['jobs']:
        label,case,method=job['label'],job['case_id'],job['method']
        result=json.loads((path/(label+'_'+method+'.json')).read_text())
        truth_path=(ROOT/'truth_development' if factor else HISTORY/'truth_after_freeze/holdout')/(case+'.npz')
        with np.load(truth_path,allow_pickle=False) as z:truth=np.asarray(z['m'],float).ravel();mask=np.asarray(z['mask'],bool).ravel()
        truths[label]=truth
        row=dict(run_name=manifest['run_name'],label=label,case_id=case,method=method,status=result['status'],
                 wall_seconds=result['wall_seconds'],cpu_seconds=result.get('cpu_seconds'),
                 scope='EXPOSED_TRUTH_DEVELOPMENT',truth_source=str(truth_path.relative_to(ROOT.parent)),
                 ranking='ONLINE_ORIGINAL_DATA_PLUS_COMMON_REGULARIZER',costcounts=result.get('costcounts'),
                 accounting=result.get('accounting'),gate=result.get('gate'),finite_endpoint=result.get('m') is not None)
        if factor:row['factor']=mapping[case]['parameters']
        if row['finite_endpoint']:
            m=np.asarray(result['m'],float).ravel();top=result['top1'];candidates=result['candidates']
            row.update(metrics(m,truth,mask),relative_data_residual=top['relative_residual'],objective=top['objective'],
                       top1_stage=top['stage'],top1_origin=top['origin'],
                       K=len(candidates),all_unknowns_released=result['all_unknowns_released'],
                       online_top1_is_score_minimum=not candidates or top['objective']<=min(c['objective'] for c in candidates)+1e-15,
                       material_min=float(m.min()),material_max=float(m.max()),
                       success_gates_pass=bool(top['relative_residual']<=.001 and np.sqrt(np.mean((m-truth)**2))<=1.5 and metrics(m,truth,mask)['IoU']>=.7))
            candmetrics=[dict(index=i,origin=c['origin'],stage=c['stage'],candidate_kind=c.get('candidate_kind'),objective=c['objective'],relative_residual=c['relative_residual'],**metrics(np.asarray(c['m'],float).ravel(),truth,mask)) for i,c in enumerate(candidates)]
            row['candidate_metrics']=candmetrics
            row['oracle_bestK_RMSE']=min((c['RMSE'] for c in candmetrics),default=row['RMSE'])
            row['oracle_interpretation']='POSTHOC_GT_RANKED; NOT ONLINE SELECTION OR CERTIFIED BRANCHES'
            row['ray_proposal_count']=len(result['raw_proposals'])
            row['ray_corrected_count']=sum(bool(q.get('corrected')) for q in result['raw_proposals'])
            row['proposal_summary']=[{key:q.get(key) for key in ['proposal_id','ray','t','raw_score','selection','corrected','corrected_score','seed']} for q in result['raw_proposals']]
            if factor:
                p=row['factor'];fg=m>=p['background']+.5*p['amplitude'];union=np.logical_or(fg,mask).sum()
                row['IoU_factor_midpoint']=float(np.logical_and(fg,mask).sum()/max(union,1))
                row['RMSE_over_amplitude']=row['RMSE']/p['amplitude']
                row['factor_success_scope']='DESCRIPTIVE_FACTOR_STUDY; fixed IoU threshold 5 is inapplicable to amplitude 1.3'
            materials[(label,method)]=m;runs[(label,method)]=result
        rows.append(row)
    scored=dict(schema='r4-development-scores-v1',run_name=manifest['run_name'],rows=rows,
                code_version_note=manifest.get('code_version_note'),
                solver_source_snapshot=manifest.get('solver_source_snapshot'),
                known_design_limits=manifest.get('known_design_limits',[]),
                score_physics_cost=dict(LU=0,complex_RHS=0),
                metric_definition=dict(IoU_foreground_threshold=5.,factor_IoU='known-factor midpoint only in posthoc factor analysis',
                    trajectory='Every completed paid full-physics trial, including rejected optimizer trials; chronological order',
                    archive='Corrected candidates sorted by original physical score; not certified branches'),
                scope='Retrospective DEVELOPMENT; historical GT was previously exposed; no fresh claim')
    dump(path/'SCORES.json',scored)
    scalar=[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in rows]
    with (path/'SCORES.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in scalar))));w.writeheader();w.writerows(scalar)
    return scored,materials,runs,truths

def figures(rows,materials,runs,truths,suffix=''):
    labels=list(dict.fromkeys(r['label'] for r in rows));figdir=ROOT/'figures';figdir.mkdir(exist_ok=True)
    methods=list(dict.fromkeys(r['method'] for r in rows))
    fig,axes=plt.subplots(len(labels),3,figsize=(14,3.6*len(labels)),squeeze=False,layout='constrained')
    for i,label in enumerate(labels):
        truth=truths[label]
        for method in methods:
            result=runs.get((label,method))
            if result is None:continue
            trace=result['trajectory'];seconds=np.asarray([q['seconds'] for q in trace]);m=np.asarray([q['m'] for q in trace]);rmse=np.sqrt(np.mean((m-truth)**2,axis=1));res=np.asarray([q['relative_residual'] for q in trace]);obj=np.asarray([q['objective'] for q in trace])
            color=COLORS[method]
            axes[i,0].semilogy(seconds,res,color=color,lw=.9,alpha=.85,label=method)
            axes[i,1].plot(seconds,rmse,color=color,lw=.9,alpha=.85)
            axes[i,2].semilogx(res,rmse,color=color,lw=.8,alpha=.7)
            axes[i,2].scatter(res[-1],rmse[-1],color=color,s=25,marker='x')
            for stage in [8,12,24]:
                indices=[k for k,q in enumerate(trace) if q['stage']==stage]
                if indices:
                    j=indices[0];axes[i,2].scatter(res[j],rmse[j],color=color,s=14,marker='o')
            # Exact plotted physical trajectories are also retained in solver JSON and NPZ.
            plotdir=ROOT/'results/plot_coordinates';plotdir.mkdir(exist_ok=True)
            np.savez_compressed(plotdir/(label+'_'+method+'_plot_coordinates.npz'),seconds=seconds,relative_residual=res,RMSE=rmse,objective=obj)
        axes[i,0].set_ylabel(label+'\nRelative data residual');axes[i,0].set_xlabel('Elapsed wall seconds')
        axes[i,1].set_ylabel('Material RMSE (posthoc GT)');axes[i,1].set_xlabel('Elapsed wall seconds')
        axes[i,2].set_ylabel('Material RMSE (posthoc GT)');axes[i,2].set_xlabel('Relative data residual')
        for ax in axes[i]:ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);fig.suptitle('R4 exposed-truth development: actual paid Maxwell trajectories\nChronological full-physics trials; circles = stage starts; x = last trial; no certified branch claim',fontsize=12)
    for extension in ['png','pdf','svg']:fig.savefig(figdir/('maxwell_branch_trajectories'+suffix+'.'+extension),dpi=180)
    plt.close(fig)
    ncols=max(3,len(methods)+1)
    fig=plt.figure(figsize=(2.5*ncols,3.1*len(labels)+3),layout='constrained')
    grid=fig.add_gridspec(len(labels)+1,ncols,height_ratios=[1.]*len(labels)+[1.1]);image=None
    for i,label in enumerate(labels):
        for j,key in enumerate(['truth']+methods):
            ax=fig.add_subplot(grid[i,j]);m=truths[label] if key=='truth' else materials.get((label,key))
            if m is not None:image=ax.imshow(m.reshape(24,24),vmin=.5,vmax=20,cmap='viridis',interpolation='nearest')
            ax.set_xticks([]);ax.set_yticks([])
            if i==0:ax.set_title('Disclosed GT' if key=='truth' else key,fontsize=10)
            if j==0:ax.set_ylabel(label)
            elif m is not None:
                r=next(r for r in rows if r['label']==label and r['method']==key)
                ax.set_xlabel('RMSE %.2f | IoU %.2f\nr %.2g | K %d'%(r['RMSE'],r['IoU'],r['relative_data_residual'],r['K']),fontsize=8)
            else:ax.set_xlabel('No finite endpoint',fontsize=8)
    boundaries=np.linspace(0,ncols,4).astype(int)
    costaxes=[fig.add_subplot(grid[-1,boundaries[j]:boundaries[j+1]]) for j in range(3)]
    names=[r['label']+'\n'+r['method'] for r in rows];colors=[COLORS.get(r['method'],'gray') for r in rows];x=np.arange(len(rows))
    values=[[r['wall_seconds'] for r in rows],[(r.get('costcounts') or {}).get('LU',0) if isinstance(r.get('costcounts'),dict) else 0 for r in rows],[(r.get('costcounts') or {}).get('actual_complex_RHS',0) if isinstance(r.get('costcounts'),dict) else 0 for r in rows]]
    for ax,v,title in zip(costaxes,values,['Inclusive solver wall seconds','Full complex LU factorizations','Actual complex right-hand sides']):
        ax.bar(x,v,color=colors);ax.set_xticks(x,names,rotation=45,ha='right',fontsize=7);ax.set_title(title,fontsize=10);ax.grid(axis='y',alpha=.2)
    fig.suptitle('R4 DEVELOPMENT: online physical-score top1 and complete solver costs\nCommon [0.5,20] scale; posthoc oracle best@K is recorded separately',fontsize=12)
    for extension in ['png','pdf','svg']:fig.savefig(figdir/('reconstruction_and_cost'+suffix+'.'+extension),dpi=180)
    plt.close(fig)

def factor_figures(rows,materials,truths):
    shapes=['DEV7','holdout7'];amplitudes=[1.3,6.,15.]
    fig=plt.figure(figsize=(14,9),layout='constrained');grid=fig.add_gridspec(3,6,height_ratios=[1,1,.95])
    colors={'DEV7':'#0072b2','holdout7':'#d55e00'}
    for i,shape in enumerate(shapes):
        for j,amplitude in enumerate(amplitudes):
            row=next(r for r in rows if r['factor']['shape']==shape and r['factor']['amplitude']==amplitude)
            label,method=row['label'],row['method']
            for offset,key in enumerate(['GT','reconstruction']):
                ax=fig.add_subplot(grid[i,j*2+offset]);m=truths[label] if key=='GT' else materials.get((label,method))
                if m is not None:ax.imshow(m.reshape(24,24),vmin=1.2,vmax=1.2+amplitude,cmap='viridis',interpolation='nearest')
                ax.set_xticks([]);ax.set_yticks([])
                ax.set_title(('Disclosed GT' if key=='GT' else 'data L-BFGS-B')+'\nAmplitude '+str(amplitude),fontsize=9)
                if j==0 and offset==0:ax.set_ylabel(shape)
                if offset:ax.set_xlabel('RMSE/A %.2f | IoU(mid) %.2f\nr %.2g'%(row['RMSE_over_amplitude'],row['IoU_factor_midpoint'],row['relative_data_residual']),fontsize=8)
    axes=[fig.add_subplot(grid[-1,0:2]),fig.add_subplot(grid[-1,2:4]),fig.add_subplot(grid[-1,4:6])]
    for shape in shapes:
        subset=sorted([r for r in rows if r['factor']['shape']==shape],key=lambda r:r['factor']['amplitude'])
        x=[r['factor']['amplitude'] for r in subset]
        for ax,key in zip(axes,['RMSE_over_amplitude','relative_data_residual','IoU_factor_midpoint']):
            ax.plot(x,[r[key] for r in subset],'-o',color=colors[shape],label=shape);ax.set_xlabel('Foreground amplitude above background 1.2');ax.set_xticks(amplitudes);ax.grid(alpha=.2)
    axes[0].set_ylabel('Material RMSE / true amplitude');axes[1].set_ylabel('Relative data residual');axes[1].set_yscale('log');axes[2].set_ylabel('IoU using true factor midpoint');axes[2].set_ylim(0,1.05);axes[0].legend(fontsize=9)
    fig.suptitle('R4 DEVELOPMENT factor separation: 2 disclosed shapes × 3 amplitudes\nSame uniform m=1 start, 20s / 900 LU / 21600 RHS cap; panel color range 1.2…1.2+A; midpoint IoU is posthoc',fontsize=12)
    for extension in ['png','pdf','svg']:fig.savefig(ROOT/'figures'/('factor_separation.'+extension),dpi=180)
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--runs',nargs='+',required=True);parser.add_argument('--factor',action='store_true');parser.add_argument('--no-figures',action='store_true');parser.add_argument('--figure-suffix',default='');args=parser.parse_args()
    rows=[];materials={};runs={};truths={};summaries=[]
    for name in args.runs:
        scored,mm,rr,tt=score_run(ROOT/'results'/name,args.factor);rows.extend(scored['rows']);materials.update(mm);runs.update(rr);truths.update(tt);summaries.append(scored)
    aggregate=ROOT/'EXPERIMENT_RESULTS.json'
    previous=json.loads(aggregate.read_text())['runs'] if aggregate.exists() else []
    runmap={r['run_name']:r for r in previous+summaries}
    allruns=list(runmap.values());allrows=[r for run in allruns for r in run['rows']]
    dump(aggregate,dict(schema='r4-experiment-results-v1',scope='EXPOSED_TRUTH_DEVELOPMENT',runs=allruns,rows=allrows,score_physics_cost=dict(LU=0,complex_RHS=0)))
    if not args.no_figures:
        if args.factor:factor_figures(rows,materials,truths)
        else:figures(rows,materials,runs,truths,args.figure_suffix)
    print(json.dumps([{k:r.get(k) for k in ['label','method','status','RMSE','IoU','relative_data_residual','K','oracle_bestK_RMSE','wall_seconds','ray_proposal_count','ray_corrected_count']} for r in rows],indent=2))

if __name__=='__main__':main()
