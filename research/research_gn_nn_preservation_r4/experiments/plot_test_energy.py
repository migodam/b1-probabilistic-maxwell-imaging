"""Plot only exported post-seal TEST energy tables. No arrays or solvers."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import argparse,csv,json,math,sys
from collections import defaultdict,Counter
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
MAIN={
    'Voxel':'v_physics_seed7_cascade2','Jfull':'j_full_seed7_cascade2',
    'Jnull':'j_null_seed7_cascade2','Gweak':'g_weak_seed7_cascade2',
    'Gdark':'g_dark_seed7_cascade2','feedback1':'recurrent_feedback1_seed7_cascade2',
    'Gweak+const':'g_weak_seed7_cascade2_const',
}
CORE=list(MAIN)[:5]
CONFIRM={name:{seed:f'{arm}_seed{seed}_cascade2'for seed in (7,17,27)}for name,arm in
         (('Voxel','v_physics'),('Jfull','j_full'),('Gweak','g_weak'))}
PHYSICAL='Maxwell_current'
AUXILIARY='auxiliary_current_head_OFFLINE_not_Maxwell_solution'
RAW='anchor_to_raw_1';FINAL='anchor_to_final_top1';AUX='anchor_to_auxiliary_head'
PRIMARY=('P_in_strong19','S_strong9','D509')
WEAK=('P_in_weak557','S_weak39','S_null528')
TITLES={'P_in_strong19':'Internal strong19','P_in_weak557':'Internal weak557',
        'S_strong9':'Receiver strong9','S_weak39':'Receiver weak39','S_null528':'Receiver null528','D509':'Dark D509'}
NOISES=(0.,.001)

def number(value):
    try:x=float(value)
    except (ValueError,TypeError):return None
    return x if math.isfinite(x)else None
def fmt(value,signed=False):
    if value is None:return 'N/R'
    return ('{:+.2e}'if signed else'{:.2e}').format(value)if value else'0'
def csvrows(path):
    if not path.is_file():return
    with path.open(newline='')as f:
        yield from csv.DictReader(f)
def mean_field(rows,field):
    if not rows:return None
    values=[number(r.get(field))for r in rows]
    return None if any(v is None for v in values)else sum(values)/len(values)

class Tables:
    def __init__(self,root,expected):
        self.root=root;self.expected=expected;self.warnings=[];self.files={}
        self.groups={};self.boot={};self.gaps=defaultdict(list);self.data=defaultdict(list);self.noise_by_obs={}
        self.notrun=[];self.receipt={};self.duplicate_boot=[]
        runs=set(MAIN.values())|{run for seeds in CONFIRM.values()for run in seeds.values()}
        required=root/'grouped_absolute_energy.csv'
        if not required.is_file():raise FileNotFoundError('NOT_READY: '+str(required)+'; no plotting or waiting was started')
        for name in ('grouped_absolute_energy.csv','mother_paired_bootstrap.csv','current_realization_gap.csv','data_residual_relation.csv','NOT_RUN.csv','RECEIPT.json'):
            self.files[name]='PRESENT'if(root/name).is_file()else'NOT_RECORDED'
        for row in csvrows(required):
            if row.get('run')not in runs:continue
            noise=number(row.get('noise'))
            if noise not in NOISES:continue
            key=(row['run'],noise,row['quantity'],row['stage'],row['projection'])
            group=self.groups.setdefault(key,dict(cases=0,sums={k:0. for k in ('before','after','delta','update','signed_cross')},missing=[]))
            count=number(row.get('cases'))
            if count is None:group['missing'].append('cases');continue
            group['cases']+=int(count)
            for field in group['sums']:
                value=number(row.get(field+'_sum'))
                if value is None:group['missing'].append(field)
                else:group['sums'][field]+=value
        for row in csvrows(root/'mother_paired_bootstrap.csv'):
            if row.get('run')not in runs or row.get('cohort')!='ALL':continue
            key=(row['run'],number(row.get('noise')),row['quantity'],row['stage'],row['projection'])
            if key in self.boot:self.duplicate_boot.append(list(key));self.boot[key]=None
            else:self.boot[key]=row
        for name,destination,stage_key in (('current_realization_gap.csv',self.gaps,False),('data_residual_relation.csv',self.data,True)):
            for row in csvrows(root/name):
                if row.get('run')not in runs:continue
                noise=number(row.get('noise'))
                if noise not in NOISES:continue
                key=(row['run'],noise,row['stage'])if stage_key else(row['run'],noise)
                destination[key].append(row);self.noise_by_obs[row['obs_id']]=noise
        self.notrun=[r for r in csvrows(root/'NOT_RUN.csv')if r.get('run')in runs]
        if(root/'RECEIPT.json').is_file():self.receipt=json.loads((root/'RECEIPT.json').read_text())
        for name in ('mother_paired_bootstrap.csv','current_realization_gap.csv','data_residual_relation.csv'):
            if self.files[name]!='PRESENT':self.warnings.append(name+': NOT_RECORDED; corresponding points/intervals remain unavailable')
        if self.duplicate_boot:self.warnings.append('Duplicate bootstrap keys are marked ambiguous; no interval chosen')
    def n(self,run,noise,stage=RAW,quantity=PHYSICAL,projection='total'):
        return self.groups.get((run,noise,quantity,stage,projection),{}).get('cases',0)
    def value(self,run,noise,projection,stage=RAW,field='after',quantity=PHYSICAL):
        group=self.groups.get((run,noise,quantity,stage,projection))
        if not group or not group['cases']or field in group['missing']or'cases'in group['missing']:return None
        return group['sums'][field]/group['cases']
    def delta(self,run,noise,projection,stage):
        row=self.boot.get((run,noise,PHYSICAL,stage,projection))
        if row is not None:
            return dict(mean=number(row.get('delta_mean')),low=number(row.get('delta_ci_low')),high=number(row.get('delta_ci_high')),
                n=number(row.get('mothers')),source='saved mother-paired bootstrap, cohort ALL')
        return dict(mean=self.value(run,noise,projection,stage,'delta'),low=None,high=None,n=self.n(run,noise,stage,projection=projection),
            source='grouped paired delta mean; bootstrap interval NOT_RECORDED')
    def count_label(self,run,noise):return f"R{self.n(run,noise,RAW)}/F{self.n(run,noise,FINAL)}"
    def realization(self,run,noise):
        gap=self.gaps.get((run,noise),[]);raw=self.data.get((run,noise,'raw_1'),[]);final=self.data.get((run,noise,'final_top1'),[])
        gids={r['obs_id']for r in gap};rids={r['obs_id']for r in raw};fids={r['obs_id']for r in final}
        duplicate=len(gids)!=len(gap)or len(rids)!=len(raw)or len(fids)!=len(final)
        matched=not duplicate and bool(gids)and gids==rids and self.n(run,noise,RAW)==len(gids)and self.n(run,noise,AUX,AUXILIARY)==len(gids)
        final_matched=matched and fids==gids and self.n(run,noise,FINAL)==len(gids)
        status='MATCHED_AUXILIARY_RAW_ROWS'if matched else'UNAVAILABLE_OR_UNMATCHED_AUXILIARY_RAW_COHORT'
        return dict(gaps=gap,raw=raw,final=final,matched=matched,final_matched=final_matched,n=len(gids),
            raw_rows=len(raw),final_rows=len(final),duplicates=duplicate,status=status)
    def counts(self):
        rows=[]
        allruns=dict(MAIN)
        for method,seeds in CONFIRM.items():
            for seed,run in seeds.items():
                if run not in allruns.values():allruns[f'{method} s{seed}']=run
        for label,run in allruns.items():
            for noise in NOISES:
                rel=self.realization(run,noise)
                nr=[r for r in self.notrun if r['run']==run and self.noise_by_obs.get(r['obs_id'])==noise]
                unknown=[r for r in self.notrun if r['run']==run and r['obs_id']not in self.noise_by_obs]
                raw=self.n(run,noise,RAW);final=self.n(run,noise,FINAL)
                rows.append(dict(method=label,run=run,noise=noise,planned_mothers=self.expected,raw_current_cases=raw,final_current_cases=final,
                    raw_coverage_gap=max(0,self.expected-raw),final_coverage_gap=max(0,self.expected-final),
                    raw_energy_value_status='RECORDED'if self.value(run,noise,'total',RAW)is not None else'NOT_RECORDED_OR_NONFINITE',
                    final_energy_value_status='RECORDED'if self.value(run,noise,'total',FINAL)is not None else'NOT_RECORDED_OR_NONFINITE',
                    auxiliary_cases=self.n(run,noise,AUX,AUXILIARY),realization_exported_rows=rel['n'],realization_matched_rows=rel['n']if rel['matched']else 0,realization_status=rel['status'],
                    raw_observed_rows=rel['raw_rows'],final_observed_rows=rel['final_rows'],final_chain_matched=rel['final_matched'],
                    auxiliary_note='N/A_NO_AUXILIARY_HEAD'if run.startswith('v_physics_')else'available only when exported',
                    NOT_RUN_events_known_noise=len(nr),NOT_RUN_reason_counts=json.dumps(Counter(r.get('reason','NOT_RECORDED')for r in nr),sort_keys=True),
                    NOT_RUN_events_noise_unassigned=len(unknown),
                    coverage='RECORDED_FULL_COUNTS'if raw==final==self.expected else'PARTIAL_OR_NOT_RECORDED',
                    interpretation='counts are availability, never scientific PASS/FAIL; all NOT_RUN events retained'))
        return rows

def axes_scale(ax,values,absolute=False):
    finite=[v for v in values if v is not None and math.isfinite(v)]
    positive=[abs(v)for v in finite if v!=0]
    if not positive:return
    if absolute and all(v>0 for v in finite):ax.set_yscale('log')
    else:ax.set_yscale('symlog',linthresh=min(positive)/5)
def missing_point(ax,x,color='0.45'):
    from matplotlib.transforms import blended_transform_factory
    ax.plot([x],[.025],marker='x',ls='',color=color,transform=blended_transform_factory(ax.transData,ax.transAxes),ms=5)
def heading(fig,subtitle,expected):
    fig.suptitle(f'NEW INDEPENDENT TEST — {expected} planned mothers per noise (clean + noisy)\n'+subtitle,fontsize=15)
def run_ticks(ax,names,runs,noise,tables):
    ax.set_xticks(range(len(names)),[name+'\n'+tables.count_label(run,noise)for name,run in zip(names,runs)],rotation=22,fontsize=8)
def save(fig,out,name,files):
    fig.savefig(out/name,dpi=170);files.append(name)
    import matplotlib.pyplot as plt
    plt.close(fig)

def single_spaces(t,out,bands,label,files):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    for delta in (False,True):
        fig,axes=plt.subplots(2,3,figsize=(19,10),layout='constrained')
        for ni,noise in enumerate(NOISES):
            for bi,band in enumerate(bands):
                ax=axes[ni,bi];values=[]
                for mi,(method,run)in enumerate(MAIN.items()):
                    stages=[(RAW,'before',-.20,'#666666','s'),(RAW,'after',0,'#e08214','^'),(FINAL,'after',.20,'#2166ac','o')]
                    if delta:stages=[(RAW,'delta',-.12,'#e08214','^'),(FINAL,'delta',.12,'#2166ac','o')]
                    for stage,field,offset,color,marker in stages:
                        estimate=t.delta(run,noise,band,stage)if delta else None
                        value=estimate['mean']if estimate else t.value(run,noise,band,stage,field)
                        if value is None:missing_point(ax,mi+offset,color);continue
                        values.append(value);ax.plot(mi+offset,value,marker=marker,color=color,ms=6,ls='')
                        if estimate and estimate['low']is not None and estimate['high']is not None:
                            ax.vlines(mi+offset,estimate['low'],estimate['high'],color=color,lw=1.2)
                            values.extend((estimate['low'],estimate['high']))
                if delta:ax.axhline(0,color='black',lw=.7)
                axes_scale(ax,values,absolute=not delta);run_ticks(ax,list(MAIN),list(MAIN.values()),noise,t)
                ax.set_title(('clean'if noise==0 else'noise .001')+' / '+TITLES[band]);ax.grid(axis='y',alpha=.2)
                ax.set_ylabel('signed mean Δ squared error vs common GN'if delta else'mean absolute squared Maxwell-current error')
        handles=[Line2D([],[],marker='s',color='#666666',ls='',label='common GN before (raw cohort)'),Line2D([],[],marker='^',color='#e08214',ls='',label='raw_1'),Line2D([],[],marker='o',color='#2166ac',ls='',label='final top1')]
        axes[0,0].legend(handles=handles[1:]if delta else handles,fontsize=8)
        heading(fig,('Saved paired changes and mother-bootstrap intervals'if delta else'Actual Maxwell-current errors GN → verified raw_1 → final')+'; fixed public τ=.1',t.expected)
        fig.get_layout_engine().set(rect=(0,.045,1,.94))
        fig.text(.5,.012,'R/F = available raw/final cases; missing values remain ×. Cohorts may differ when failures occur. Discrete squared norms, no energy shares. '+('Intervals are per seed and available paired mothers; no new bootstrap.'if delta else'GN points use the raw cohort; final has its displayed own coverage.'),ha='center',fontsize=9)
        save(fig,out,label+('_paired_delta.png'if delta else'_absolute.png'),files)

def twofold(t,out,files):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    definitions=('intersection','I_to_O','O_to_I')
    titles=('Numerical intersections + remainder47','Sequential I→O: P_O P_I','Sequential O→I: P_I P_O')
    for noise in NOISES:
        valid=[abs(t.value(run,noise,f'{d}__{i}__{o}',FINAL,'delta')or 0)for run in MAIN.values()for d in definitions for i in ('strong','weak')for o in ('strong','weak','null')]
        positive=[v for v in valid if v>0];vmax=max(positive)if positive else 1.;lin=min(positive)/5 if positive else .1
        norm=SymLogNorm(linthresh=lin,vmin=-vmax,vmax=vmax)
        fig,axes=plt.subplots(len(MAIN),3,figsize=(20,26));fig.subplots_adjust(left=.12,right=.98,bottom=.067,top=.91,hspace=.78,wspace=.30)
        image=None
        for ri,(method,run)in enumerate(MAIN.items()):
            for ci,definition in enumerate(definitions):
                ax=axes[ri,ci];color=np.full((2,3),np.nan)
                for i,ib in enumerate(('strong','weak')):
                    for j,ob in enumerate(('strong','weak','null')):
                        value=t.value(run,noise,f'{definition}__{ib}__{ob}',FINAL,'delta')
                        if value is not None:color[i,j]=value
                image=ax.imshow(np.ma.masked_invalid(color),norm=norm,cmap='RdBu_r',aspect='auto')
                ax.set_facecolor('.88');ax.set_xticks(range(3),('S strong9','S weak39','S null528'),fontsize=9);ax.xaxis.tick_top()
                ax.set_yticks(range(2),('I strong19','I weak557'),fontsize=9)
                if ri==0:ax.set_title(titles[ci],pad=30,fontsize=12)
                for i,ib in enumerate(('strong','weak')):
                    for j,ob in enumerate(('strong','weak','null')):
                        band=f'{definition}__{ib}__{ob}'
                        text='GN_R '+fmt(t.value(run,noise,band,RAW,'before'))+'\nGN_F '+fmt(t.value(run,noise,band,FINAL,'before'))+'\nΔraw '+fmt(t.value(run,noise,band,RAW,'delta'),True)+'\nΔfinal '+fmt(t.value(run,noise,band,FINAL,'delta'),True)
                        bg=image.cmap(image.norm(color[i,j]))if math.isfinite(color[i,j])else(1,1,1,1)
                        lum=.2126*bg[0]+.7152*bg[1]+.0722*bg[2]
                        ax.text(j,i,text,ha='center',va='center',fontsize=7.8,color='white'if lum<.45 else'black')
                if ci==0:
                    ax.text(-.34,.5,method+'\n'+t.count_label(run,noise),transform=ax.transAxes,ha='center',va='center',fontsize=10,fontweight='bold')
                    vals=[t.value(run,noise,'intersection_remainder',RAW,'before'),t.value(run,noise,'intersection_remainder',FINAL,'before'),t.value(run,noise,'intersection_remainder',RAW,'delta'),t.value(run,noise,'intersection_remainder',FINAL,'delta')]
                    ax.text(.5,-.10,'Remainder47: GN_R '+fmt(vals[0])+'; GN_F '+fmt(vals[1])+'\nΔraw '+fmt(vals[2],True)+'; Δfinal '+fmt(vals[3],True),transform=ax.transAxes,ha='center',va='top',fontsize=8)
                else:ax.text(.5,-.10,'Sequential energies; not orthogonal joint intersection cells',transform=ax.transAxes,ha='center',va='top',fontsize=8)
        heading(fig,('clean'if noise==0 else'noise .001')+'; different Twofold definitions, actual Maxwell-current errors; public τ=.1',t.expected)
        fig.text(.5,.941,'Cell text retains both original GN cohort means and signed raw/final changes. No forced 100% allocation or share normalization.',ha='center',fontsize=11)
        cax=fig.add_axes((.25,.035,.50,.010));fig.colorbar(image,cax=cax,orientation='horizontal',label='signed Δfinal; symmetric-log color display only')
        fig.text(.5,.009,'Numerical intersection ranks [0,0,0]/[0,20,509], with orthogonal remainder47 shown separately. R/F coverage and NOT_RUN are retained.',ha='center',fontsize=10)
        save(fig,out,'twofold_'+('clean'if noise==0 else'noisy')+'.png',files)

def realization(t,out,files):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    colors=dict(zip(MAIN,('#777777','#377eb8','#ff7f00','#4daf4a','#984ea3','#e41a1c','#a65628')))
    bands=('total',)+PRIMARY
    fig,axes=plt.subplots(2,4,figsize=(22,10),layout='constrained');unavailable=[]
    for ni,noise in enumerate(NOISES):
        for bi,band in enumerate(bands):
            ax=axes[ni,bi];values=[]
            for method,run in MAIN.items():
                rel=t.realization(run,noise)
                if not rel['matched']:
                    if bi==0:unavailable.append(dict(method=method,noise=noise,status='N/A_NO_AUXILIARY_HEAD'if method=='Voxel'else rel['status']))
                    continue
                if bi==0 and not rel['final_matched']:unavailable.append(dict(method=method,noise=noise,status='FINAL_COHORT_UNMATCHED_POINT_NOT_CONNECTED'))
                stages=[t.value(run,noise,band,RAW,'before'),t.value(run,noise,band,AUX,'after',AUXILIARY),t.value(run,noise,band,RAW,'after'),t.value(run,noise,band,FINAL,'after')if rel['final_matched']else None]
                y=[float('nan')if v is None else v for v in stages];values.extend(v for v in stages if v is not None)
                ax.plot(range(4),y,'o-',color=colors[method],label=method+f" n={rel['n']}",markersize=4,linewidth=1.2)
            axes_scale(ax,values,absolute=True);ax.set_xticks(range(4),('GN\nphysical','auxiliary\nhead','raw_1\nphysical','final\nphysical'),fontsize=8)
            ax.set_title(('clean'if noise==0 else'noise .001')+' / '+TITLES.get(band,'Total current discrepancy'));ax.grid(alpha=.2)
            if bi==0:
                ax.set_ylabel('mean discrete squared current discrepancy')
                if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=7)
    heading(fig,'Matched NN auxiliary head → physical realization → final; public τ=.1',t.expected)
    fig.get_layout_engine().set(rect=(0,.10,1,.93))
    missing_text='; '.join(('clean 'if r['noise']==0 else'noisy ')+r['method']+(' final'if r['status'].startswith('FINAL')else' aux/raw')for r in unavailable)or'none'
    fig.text(.5,.070,'Unavailable/unmatched, explicitly retained: '+missing_text,ha='center',fontsize=8)
    fig.text(.5,.035,'Auxiliary points are head discrepancies, never called true Maxwell current errors. Voxel has no auxiliary head (N/A).',ha='center',fontsize=10)
    fig.text(.5,.012,'Only matched NN auxiliary/raw cohorts are connected. Final is omitted when its cohort differs; no baseline incumbent is assumed to share an auxiliary state.',ha='center',fontsize=9)
    save(fig,out,'auxiliary_to_physical_chain.png',files)
    nn=[(name,run)for name,run in MAIN.items()if name!='Voxel']
    fig,axes=plt.subplots(2,3,figsize=(24,11),layout='constrained')
    fieldsets=[('State / realization gap (not GT error)',('state_defect_energy','current_realization_gap_energy','internal_realization_gap_energy'),('state defect','current gap','internal-field gap')),
               ('Receiver update relative to common GN',('aux_receiver_update_energy','receiver_realization_gap_energy','signed_receiver_cross','physical_receiver_update_energy'),('aux update','realization gap','signed cross','physical update')),
               ('Raw physical receiver error / misfit / noise',('true_S_field_error_energy','observed_misfit_energy','noise_energy','signed_residual_noise_cross'),('true error','observed misfit','noise','signed cross'))]
    for ni,noise in enumerate(NOISES):
        for ci,(title,fields,labels)in enumerate(fieldsets):
            ax=axes[ni,ci];values=[]
            for mi,(method,run)in enumerate(nn):
                rel=t.realization(run,noise)
                rows=rel['gaps']if ci<2 else rel['raw']
                for fi,(field,label)in enumerate(zip(fields,labels)):
                    value=mean_field(rows,field)if rel['matched']else None
                    x=mi+(fi-(len(fields)-1)/2)*.17
                    if value is None:missing_point(ax,x);continue
                    ax.bar(x,value,width=.16,color=f'C{fi}',label=label if mi==0 else None);values.append(value)
            ax.axhline(0,color='black',lw=.7);axes_scale(ax,values,absolute=ci==0)
            ax.set_xticks(range(len(nn)),[name+'\n'+('paired n='+str(t.realization(run,noise)['n'])if t.realization(run,noise)['matched']else'N/R; rows='+str(t.realization(run,noise)['n']))for name,run in nn],rotation=22,fontsize=8)
            ax.set_title(('clean'if noise==0 else'noise .001')+' / '+title);ax.grid(axis='y',alpha=.2)
            ax.set_ylabel('mean discrete squared norm / signed cross');ax.legend(handles=[Patch(color=f'C{fi}',label=label)for fi,label in enumerate(labels)],fontsize=8)
    heading(fig,'NN matched raw-state realization and receiver accounting; saved scalars only',t.expected)
    fig.get_layout_engine().set(rect=(0,.065,1,.93))
    fig.text(.5,.035,'Negative cross terms remain signed. These are separate saved quantities; no closure/identity check, re-solve, or 100% allocation is performed.',ha='center',fontsize=10)
    fig.text(.5,.012,'Receiver-update panel is relative to GN, not GT; receiver-error panel compares raw physical output with GT / observed Y as labelled. Missing/unmatched cohorts remain ×.',ha='center',fontsize=9)
    save(fig,out,'realization_state_receiver_accounting.png',files)
    return unavailable

def confirmation(t,out,files):
    import matplotlib.pyplot as plt
    colors={7:'#377eb8',17:'#ff7f00',27:'#4daf4a'}
    combinations=[(method,seed,run)for method,seeds in CONFIRM.items()for seed,run in seeds.items()]
    fig,axes=plt.subplots(3,2,figsize=(19,14),layout='constrained')
    for bi,band in enumerate(PRIMARY):
        for ni,noise in enumerate(NOISES):
            ax=axes[bi,ni];values=[]
            for mi,(method,seed,run)in enumerate(combinations):
                for stage,offset,marker in ((RAW,-.12,'^'),(FINAL,.12,'o')):
                    e=t.delta(run,noise,band,stage);v=e['mean']
                    if v is None:missing_point(ax,mi+offset,colors[seed]);continue
                    ax.plot(mi+offset,v,marker=marker,color=colors[seed],ls='',mfc='white'if stage==RAW else colors[seed]);values.append(v)
                    if e['low']is not None and e['high']is not None:ax.vlines(mi+offset,e['low'],e['high'],color=colors[seed]);values.extend((e['low'],e['high']))
            ax.axhline(0,color='black',lw=.7);axes_scale(ax,values)
            ax.set_xticks(range(len(combinations)),[f'{m} s{s}\n'+t.count_label(run,noise)for m,s,run in combinations],rotation=25,fontsize=8)
            ax.set_title(('clean'if noise==0 else'noise .001')+' / '+TITLES[band]);ax.set_ylabel('paired signed Δ error versus common GN');ax.grid(axis='y',alpha=.2)
    heading(fig,'Three-seed confirmation kept separate: Voxel, Jfull, Gweak; public τ=.1',t.expected)
    fig.get_layout_engine().set(rect=(0,.06,1,.94))
    fig.text(.5,.025,'Open triangles: raw; filled circles: final. Each interval uses its own seed’s paired mothers. No seed pooling, seed-as-mother counting, or pooled seed standard error.',ha='center',fontsize=10)
    save(fig,out,'three_seed_confirmation.png',files)

def sensitivity(t,out,files):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(19,10),layout='constrained')
    categories=(('Gin_strong','P_in_strong19'),('S_strong','S_strong9'),('dark','D509'))
    for ni,noise in enumerate(NOISES):
        for bi,(suffix,base)in enumerate(categories):
            ax=axes[ni,bi];values=[]
            for mi,name in enumerate(CORE):
                for ti,tau in enumerate((.03,.1,.3)):
                    projection=base if tau==.1 else f'tau{tau:g}__{suffix}'
                    value=t.value(MAIN[name],noise,projection,FINAL,'delta')
                    if value is None:missing_point(ax,mi+(ti-1)*.16,f'C{ti}');continue
                    ax.plot(mi+(ti-1)*.16,value,'o',color=f'C{ti}',label=f'τ={tau}'if mi==0 else None);values.append(value)
            ax.axhline(0,color='black',lw=.7);axes_scale(ax,values)
            ax.set_xticks(range(len(CORE)),[name+'\n'+t.count_label(MAIN[name],noise)for name in CORE],rotation=20,fontsize=8)
            ax.set_title(('clean'if noise==0 else'noise .001')+' / '+suffix);ax.set_ylabel('paired mean Δfinal in the stated public space');ax.grid(alpha=.2)
            from matplotlib.lines import Line2D
            ax.legend(handles=[Line2D([],[],marker='o',color=f'C{ti}',ls='',label=f'τ={tau}')for ti,tau in enumerate((.03,.1,.3))],fontsize=8)
    heading(fig,'Core public-threshold sensitivity; saved descriptive means, no new fitting',t.expected)
    fig.get_layout_engine().set(rect=(0,.065,1,.94))
    fig.text(.5,.035,'Different τ values define different ranks/spaces. No sensitivity bootstrap intervals were exported; no τ is selected by these plots.',ha='center',fontsize=10)
    fig.text(.5,.012,'Training and methods remain frozen. τ=.1 is primary; .03/.3 are post-seal projection sensitivity, not additional trained models.',ha='center',fontsize=9)
    save(fig,out,'core_threshold_sensitivity.png',files)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=HERE/'results/test_energy')
    p.add_argument('--output-dir',type=Path,default=HERE/'figures/test_energy')
    p.add_argument('--no-sensitivity',action='store_true')
    a=p.parse_args()
    try:t=Tables(a.input_dir,108)
    except FileNotFoundError as error:print(str(error),file=sys.stderr);return 2
    # No score, GT array, operator, model, NPZ or solver is imported/read.
    import matplotlib
    matplotlib.use('Agg')
    a.output_dir.mkdir(parents=True,exist_ok=True);files=[]
    single_spaces(t,a.output_dir,PRIMARY,'strong_dark',files)
    single_spaces(t,a.output_dir,WEAK,'weak_null',files)
    twofold(t,a.output_dir,files);unavailable=realization(t,a.output_dir,files)
    confirmation(t,a.output_dir,files)
    if not a.no_sensitivity:sensitivity(t,a.output_dir,files)
    counts=t.counts()
    with(a.output_dir/'method_counts.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(counts[0]));w.writeheader();w.writerows(counts)
    primary=[]
    for method,run in MAIN.items():
        for noise in NOISES:
            for band in PRIMARY+WEAK:
                primary.append(dict(method=method,run=run,noise=noise,projection=band,
                    GN_raw_cohort=t.value(run,noise,band,RAW,'before'),raw=t.value(run,noise,band,RAW,'after'),
                    GN_final_cohort=t.value(run,noise,band,FINAL,'before'),final=t.value(run,noise,band,FINAL,'after'),
                    raw_delta=t.delta(run,noise,band,RAW),final_delta=t.delta(run,noise,band,FINAL)))
    summary=dict(status='PLOTTED_EXPORTED_NEW_INDEPENDENT_TEST_TABLES',expected_mothers_per_noise=108,
        source_files=t.files,producer_receipt_status=t.receipt.get('status','NOT_RECORDED'),
        before_reference=t.receipt.get('before_reference','common strict GN per exporting code; producer metadata NOT_RECORDED'),
        main_registry=MAIN,confirmation_registry=CONFIRM,method_counts=counts,primary_energies=primary,
        NOT_RUN_events=len(t.notrun),NOT_RUN_reasons=dict(Counter(r.get('reason','NOT_RECORDED')for r in t.notrun)),
        groups_with_missing_or_nonfinite_saved_fields=sum(bool(g['missing'])for g in t.groups.values()),
        missing_or_nonfinite_field_counts=dict(Counter(field for group in t.groups.values()for field in group['missing'])),
        auxiliary_chain_unavailable=unavailable,warnings=t.warnings,duplicate_bootstrap_keys=t.duplicate_boot,
        figures=files,geometry_metadata=t.receipt.get('geometry','NOT_RECORDED'),
        interpretation='Availability counts are not scientific PASS/FAIL. Actual final subsets and NOT_RUN remain explicit; never infer complete coverage from the title.',
        statistical_scope='source mother-paired bootstrap within one seed and noise, cohort ALL; no pooling seeds as independent mothers; no new bootstrap',
        twofold_scope='genuine numerical intersections plus remainder47 are separate from both sequential orders; no forced 100% allocation',
        auxiliary_scope='only matched known-NN auxiliary/raw rows, never baseline-incumbent pairing; auxiliary discrepancies are not true Maxwell-solution current errors',
        sensitivity_scope='tau .1 primary, .03/.3 descriptive saved projection means without exported CI; different ranks, no method or threshold selection',
        scientific_judgment='parent interpretation after finite TEST seal; no plotting acceptance thresholds or claims of superiority',
        physics_calls=0,GT_array_reads=0,TEST_score_reads=0,formula_identity_checks=0)
    (a.output_dir/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=summary['status'],figures=len(files),method_count_rows=len(counts),output_dir=str(a.output_dir),warnings=t.warnings),indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
