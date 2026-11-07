"""Two scalar-only ablation figures; no GT arrays, physics or bootstrap."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import json,math,statistics
from collections import defaultdict,Counter
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parents[1]
OUT=HERE/'figures/ablations'
ALPHAS={0.:'j_full_seed7_cascade2_posthoc_a0',.1:'j_full_seed7_cascade2_posthoc_a0.1',
        .3:'j_full_seed7_cascade2_posthoc_a0.3',.6:'j_full_seed7_cascade2_posthoc_a0.6',1.:'j_full_seed7_cascade2'}
PANELS=[('Jfull loss ablation',[
    ('material','j_full_seed7_cascade2'),('material +\nconstitutive','j_full_seed7_cascade2_const'),('material +\ninternal','j_full_seed7_cascade2_internal')]),
    ('Gweak loss ablation',[
    ('material','g_weak_seed7_cascade2'),('material +\nconstitutive','g_weak_seed7_cascade2_const'),('material +\ninternal','g_weak_seed7_cascade2_internal')]),
    ('Matched FT from selected Gweak+const',[
    ('selected const\nno FT','g_weak_seed7_cascade2_const'),('material FT\n+320 updates','g_weak_seed7_cascade2_const_ft_material'),('full-VIE FT\n+320 updates','g_weak_seed7_cascade2_const_ft_fullvie')])]

def finite(x):return isinstance(x,(int,float))and math.isfinite(x)
def val_data(scores):
    result=[]
    for noise in (0.,.001):
        for alpha,run in ALPHAS.items():
            rows=[r for r in scores['rows']if r.get('split')=='val'and r.get('run')==run and float(r['noise'])==noise]
            parents=defaultdict(list)
            for row in rows:
                top=row.get('top1');value=top.get('NRMSE')if isinstance(top,dict)else None
                if finite(value):parents[row['parent_id']].append(float(value))
            means={parent:statistics.fmean(values)for parent,values in parents.items()}
            planned={r['parent_id']for r in rows}
            qualified=sum(r.get('selection_status')=='PASS'and isinstance(r.get('top1'),dict)and r['top1'].get('physical_data_feasible',False)for r in rows)
            result.append(dict(alpha=alpha,run=run,noise=noise,planned_mothers=78,recorded_mothers=len(planned),
                parents_with_NRMSE=len(means),parent_mean_NRMSE=statistics.fmean(means.values())if means else None,
                physically_qualified_count=qualified,joint_success_count=sum(r.get('joint_success')is True for r in rows),
                joint_success_total55_count=sum(r.get('joint_success_total55')is True for r in rows),
                missing_top1_count=sum(not isinstance(r.get('top1'),dict)for r in rows),
                status_counts=dict(Counter(r.get('status','NOT_RECORDED')for r in rows)),
                anchor_fallback_count=sum(isinstance(r.get('top1'),dict)and r['top1'].get('label')=='anchor'for r in rows),
                selected_endpoint=alpha==1.,endpoint1_reuses_original_Jfull=alpha==1.))
    return result

def plot_val(rows):
    fig,axes=plt.subplots(2,2,figsize=(15,10))
    fig.subplots_adjust(left=.08,right=.97,bottom=.14,top=.84,hspace=.37,wspace=.25)
    for ni,noise in enumerate((0.,.001)):
        data=[r for r in rows if r['noise']==noise];x=[r['alpha']for r in data]
        ax=axes[ni,0];y=[r['parent_mean_NRMSE']for r in data]
        ax.plot(x,y,'o-',color='#2166ac',lw=1.8)
        for xx,yy in zip(x,y):
            if yy is not None:ax.annotate(f'{yy:.5f}',(xx,yy),xytext=(0,-18)if xx==0 else(0,9),textcoords='offset points',ha='center',fontsize=9)
        last=data[-1]
        if last['parent_mean_NRMSE']is not None:ax.plot([1.],[last['parent_mean_NRMSE']],'*',color='#762a83',ms=16,zorder=4)
        available=[v for v in y if v is not None]
        if available:
            span=max(available)-min(available);pad=max(span*.18,max(available)*.035)
            ax.set_ylim(max(0,min(available)-pad),max(available)+pad*1.65)
        ax.set_ylabel('mean selected-top1 material NRMSE\n(equal weight per parent)')
        ax.set_title(('clean'if noise==0 else'noise .001')+' / material quality')
        ax=axes[ni,1]
        ax.plot(x,[r['physically_qualified_count']for r in data],'s--',color='#4d4d4d',label='physical qualification')
        ax.plot(x,[r['joint_success_count']for r in data],'o-',color='#1b7837',label='joint success')
        for xx,r in zip(x,data):
            ax.annotate(str(r['physically_qualified_count'])+'/78',(xx,r['physically_qualified_count']),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8,color='#4d4d4d')
            ax.annotate(str(r['joint_success_count'])+'/78',(xx,r['joint_success_count']),xytext=(0,8),textcoords='offset points',ha='center',fontsize=9,color='#1b7837')
        ax.set_ylim(-4,88);ax.set_ylabel('individual parent count / 78');ax.set_title(('clean'if noise==0 else'noise .001')+' / qualification');ax.legend(fontsize=9,loc='center right')
        for ax in axes[ni]:
            ax.axvline(1.,color='#762a83',ls=':',lw=1)
            ax.set_xticks(list(ALPHAS),('0','.1','.3','.6','1*'))
            ax.set_xlim(-.045,1.07);ax.set_xlabel('posthoc projection α');ax.grid(alpha=.2)
    fig.suptitle('VAL POSTHOC PROJECTION — NOT trained SoftDark\nFive fixed α values; selected endpoint α=1 reuses original Jfull',fontsize=16,y=.965)
    fig.text(.5,.885,'Seed7, two passes, 78 mothers per noise. All five points retained; no new training and no bootstrap.',ha='center',fontsize=11)
    fig.text(.5,.074,'* α=1 is the existing Jfull endpoint, not a separately trained SoftDark model. Physical qualification and joint success are distinct counts.',ha='center',fontsize=10)
    fig.text(.5,.038,'NRMSE is the available-parent mean; it is not a mean-gate success decision. Counts use all 78 planned parents, including any failures.',ha='center',fontsize=10)
    fig.savefig(OUT/'val_posthoc_alpha.png',dpi=180);plt.close(fig)

def test_data(summary):
    wanted={run for title,panel in PANELS for label,run in panel}
    rows=[r for r in summary['results']if r.get('run')in wanted and r.get('group')=='ALL']
    indexed={(r['run'],float(r['noise'])):r for r in rows}
    output=[]
    for noise in (0.,.001):
        for title,panel in PANELS:
            for label,run in panel:
                source=indexed.get((run,noise),{})
                output.append(dict(panel=title,label=label,run=run,noise=noise,
                    planned_mothers=source.get('planned_observations',108),independent_mothers=source.get('independent_mothers'),
                    saved=source.get('saved',0),qualified=source.get('qualified',0),not_run=source.get('not_run'),
                    mean_NRMSE=source.get('top1_NRMSE_mean'),joint_success=source.get('joint_success',0),
                    joint_success_total55=source.get('joint_success_total55'),method_status_counts=source.get('method_status_counts',{}),
                    source_status='RECORDED'if source else'NOT_RECORDED',
                    interpretation='mean NRMSE and per-parent joint success are separate outcomes'))
    return output

def plot_test(rows):
    fig,axes=plt.subplots(2,3,figsize=(19,11))
    fig.subplots_adjust(left=.065,right=.98,bottom=.20,top=.84,hspace=.43,wspace=.28)
    colors=('#377eb8','#e08214','#4daf4a')
    for ni,noise in enumerate((0.,.001)):
        for pi,(title,panel)in enumerate(PANELS):
            ax=axes[ni,pi];selected=[r for r in rows if r['noise']==noise and r['panel']==title]
            vals=[r['mean_NRMSE']for r in selected];finite_vals=[v for v in vals if finite(v)]
            upper=max(finite_vals)*1.36 if finite_vals else 1.
            for i,r in enumerate(selected):
                value=r['mean_NRMSE']
                if finite(value):
                    ax.bar(i,value,color=colors[i],width=.62)
                    ax.text(i,value+upper*.025,f"NRMSE {value:.5f}\njoint {r['joint_success']}/{r['planned_mothers']}",ha='center',va='bottom',fontsize=10)
                    ax.text(i,upper*.018,f"qualified {r['qualified']}/{r['planned_mothers']}",ha='center',va='bottom',fontsize=8,color='white'if value>upper*.13 else'black')
                else:ax.text(i,upper*.2,'NOT_RECORDED',ha='center',rotation=90,fontsize=9)
            ax.set_ylim(0,upper);ax.set_xticks(range(3),[r['label']for r in selected],fontsize=10)
            ax.set_title(('clean'if noise==0 else'noise .001')+' / '+title,fontsize=12)
            ax.set_ylabel('mean selected-top1 material NRMSE');ax.grid(axis='y',alpha=.2)
    fig.suptitle('NEW INDEPENDENT TEST — loss and matched-fine-tuning ablations\nSeed7, two passes; all 108 mothers per noise; every requested branch retained',fontsize=16,y=.967)
    fig.text(.5,.885,'Constitutive / internal loss weight λ=0.01. Joint labels are individual success counts, not decisions inferred from a mean.',ha='center',fontsize=11)
    fig.text(.5,.135,'Matched FT starts from the same selected Gweak+const checkpoint: 320 additional updates per FT arm.',ha='center',fontsize=11)
    fig.text(.5,.096,'Recorded additional OFFLINE FT work: 2,560 LU / 61,440 complex RHS. Training cost is separate from online inference/postprocessing budgets.',ha='center',fontsize=11)
    fig.text(.5,.057,'Mean NRMSE, physical qualification and joint success remain distinct. No mean-threshold gate, new bootstrap, or omission of worse branches.',ha='center',fontsize=10)
    fig.savefig(OUT/'test_loss_quality.png',dpi=180);plt.close(fig)

def main():
    val_path=HERE/'results/VAL_ALPHA_SCORES.json';test_path=HERE/'results/RESULTS_SUMMARY.json'
    val=json.loads(val_path.read_text());test=json.loads(test_path.read_text())
    vr=val_data(val);tr=test_data(test);OUT.mkdir(parents=True,exist_ok=True)
    plot_val(vr);plot_test(tr)
    common=dict(physics_calls=0,raw_GT_reads=0,training_calls=0,new_bootstrap=False,CPU_workers=1,BLAS_threads=1)
    (OUT/'val_posthoc_alpha.json').write_text(json.dumps(dict(common,scope='VAL_MODEL_SELECTION_POSTHOC_PROJECTION_NOT_TRAINED_SOFTDARK',source=str(val_path),alphas=list(ALPHAS),selected_alpha=1.,endpoint1='reuse original j_full_seed7_cascade2; no additional trained weight set',data=vr),indent=2,allow_nan=False)+'\n')
    (OUT/'test_loss_quality.json').write_text(json.dumps(dict(common,scope='NEW_INDEPENDENT_TEST_SCALAR_LOSS_ABLATION',source=str(test_path),source_status=test.get('status'),
        expected_mothers_per_noise=108,loss_strength=.01,FT_updates_per_arm=320,
        recorded_additional_offline_FT_work=dict(LU=2560,complex_RHS=61440,source='parent supplied actual cost in task; compact score summary has no per-arm allocation',included_in_online_budget=False),
        no_mean_gate_inference=True,all_requested_branches_retained=True,data=tr),indent=2,allow_nan=False)+'\n')
    print('\n'.join(str(OUT/name)for name in ('val_posthoc_alpha.png','val_posthoc_alpha.json','test_loss_quality.png','test_loss_quality.json')))

if __name__=='__main__':main()
