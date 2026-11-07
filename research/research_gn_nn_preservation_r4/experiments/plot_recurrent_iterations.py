"""Plot saved recurrent-depth VAL CSVs and receipts only; never run physics."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import argparse,csv,json,math,sys
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
DEFAULT_RUNS={
    'refreshed_q':'recurrent_refreshed_q_seed7_cascade2',
    'feedback1':'recurrent_feedback1_seed7_cascade2',
    'feedback_hinge':'recurrent_feedback_hinge_seed7_cascade2',
}
LABELS={'refreshed_q':'refreshed-q','feedback1':'quadratic feedback','feedback_hinge':'hinge feedback'}
COLORS={'refreshed_q':'#756bb1','feedback1':'#d95f02','feedback_hinge':'#1b9e77'}
METRICS={
    'material_RMS_step':('Material update RMS','contrast RMS'),
    'current_relative_step':('Auxiliary-current update','relative to ‖j₀‖'),
    'constitutive_relative_defect':('Constitutive defect','relative to ‖j₀‖'),
    'true_FullVIE_relative_data':('Actual Maxwell data residual','relative to ‖Y‖'),
    'physical_current_gap_relative':('Physical / auxiliary current gap','relative to ‖j₀‖'),
    'physical_internal_field_gap_relative':('Physical / network internal-field gap','relative to ‖Q₀‖'),
}

def parsed(value):
    if value is None or str(value).strip()=='':return None,'NOT_RECORDED'
    try:x=float(value)
    except (TypeError,ValueError):return None,'NONNUMERIC_RECORDED'
    return (x,'RECORDED')if math.isfinite(x)else(None,'NONFINITE_RECORDED')

def input_plan(root,mapping,contract):
    # Inventory the complete finite plan before any output is written.
    paths={mode:root/directory for mode,directory in mapping.items()}
    absent=[str(base/name)for base in paths.values()for name in ('trajectory.csv','RECEIPT.json')if not(base/name).is_file()]
    if absent:
        raise FileNotFoundError('NOT_READY: saved diagnostic inputs are absent; no waiting or numerical launch is performed.\n'+'\n'.join(absent))
    expected={c['obs_id']:c for c in contract['cases']}
    if len(expected)!=6 or any(c.get('split')!='val'for c in expected.values()):raise ValueError('The predeclared six-case VAL contract is required')
    datasets={};receipts={}
    for mode,base in paths.items():
        receipt=json.loads((base/'RECEIPT.json').read_text())
        if receipt.get('scope')!='VAL_ONLY_PUBLIC_RESIDUALS_NO_GT_SCORING':raise ValueError('Receipt is not explicitly public-input VAL: '+str(base))
        if int(receipt['trained_iterations'])!=int(contract['trained_iterations']):raise ValueError('Trained-depth metadata changed')
        if int(receipt['diagnostic_iterations'])!=int(contract['iterations']):raise ValueError('Diagnostic-depth metadata changed')
        table={}
        with(base/'trajectory.csv').open()as f:
            reader=csv.DictReader(f)
            missing_fields=set(METRICS)-set(reader.fieldnames or[])
            if missing_fields:raise ValueError('Missing saved metric columns: '+str(sorted(missing_fields)))
            for row in reader:
                if row.get('split')!='val':raise ValueError('Non-VAL row in saved trajectory')
                obs=row['obs_id'];iteration=int(row['iteration'])
                if obs not in expected:raise ValueError('Unplanned observation in saved trajectory: '+obs)
                if row['parent_id']!=expected[obs]['parent_id']or float(row['noise'])!=float(expected[obs]['noise']):raise ValueError('Trajectory identity differs from predeclared case')
                if not 1<=iteration<=int(contract['iterations']):raise ValueError('Out-of-plan recorded iteration')
                if (obs,iteration)in table:raise ValueError('Duplicate observation/iteration row')
                table[obs,iteration]=row
        datasets[mode]=table;receipts[mode]=receipt
    return paths,datasets,receipts

def endpoints(mapping,paths,datasets,receipts,contract):
    first=int(contract['trained_iterations']);last=int(contract['iterations']);flat=[];cases=[]
    for mode in mapping:
        table=datasets[mode]
        for case in contract['cases']:
            obs=case['obs_id'];record={**case,'mode':mode,'run_directory':str(paths[mode]),
                'trained_depth':first,'extended_depth':last,'receipt_status':receipts[mode].get('status','NOT_RECORDED'),
                'recorded_iterations':[i for i in range(1,last+1)if(obs,i)in table],
                'missing_iterations':[i for i in range(1,last+1)if(obs,i)not in table],'metrics':{}}
            for metric in METRICS:
                left=table.get((obs,first),{}).get(metric);right=table.get((obs,last),{}).get(metric)
                a,astatus=parsed(left);b,bstatus=parsed(right)
                delta=b-a if a is not None and b is not None else None
                record['metrics'][metric]=dict(iteration2=a,iteration16=b,delta16_minus2=delta,
                    iteration2_status=astatus,iteration16_status=bstatus,iteration2_raw=left,iteration16_raw=right)
                flat.append(dict(mode=mode,run_directory=str(paths[mode]),obs_id=obs,parent_id=case['parent_id'],family=case['family'],digit=case.get('digit'),noise=case['noise'],amplitude=case['amplitude'],
                    metric=metric,iteration2=a,iteration16=b,delta16_minus2=delta,
                    iteration2_status=astatus,iteration16_status=bstatus,iteration2_raw=left,iteration16_raw=right,
                    missing_iterations=json.dumps(record['missing_iterations']),receipt_status=record['receipt_status']))
            cases.append(record)
    provenance={mode:dict(trajectory=str(paths[mode]/'trajectory.csv'),receipt=str(paths[mode]/'RECEIPT.json'),
        checkpoint=receipts[mode].get('checkpoint'),status=receipts[mode].get('status'),
        wall_seconds=receipts[mode].get('wall_seconds'),NN_inference_seconds=receipts[mode].get('NN_inference_seconds'),
        physical_costcounts=receipts[mode].get('physical_costcounts'),
        mode_label_source='explicit directory mapping; no checkpoint is loaded by plotting')for mode in mapping}
    document=dict(status='SAVED_VAL_DEPTH_DIAGNOSTIC_PLOTTED',scope='VAL_PUBLIC_INPUT_OUT_OF_TRAINING_DEPTH_NOT_TEST',
        trained_depth=first,extended_depth=last,predeclared_cases=contract['cases'],
        selection=contract.get('selection'),inputs=provenance,endpoints=cases,
        interpretation='2-to-16 endpoint changes and finite trajectories are descriptive; no fixed-point, GT-convergence, contraction or stability guarantee',
        numerical_calls_by_this_plotter=0,GT_files_read=0,training_calls=0,
        endpoint_statuses='RECORDED/NOT_RECORDED/NONFINITE_RECORDED are data availability labels, not scientific PASS/FAIL',
        normalization_note='constitutive and current-gap denominators are saved anchor ||j0||, field-gap denominator ||Q0||; true_FullVIE_relative_data is actual-model residual against observed Y, not GT error')
    return flat,document

def chart(mapping,datasets,receipts,contract,output):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.transforms import blended_transform_factory
    nrows=len(contract['cases']);ncols=len(METRICS);last=int(contract['iterations']);trained=int(contract['trained_iterations'])
    fig,axes=plt.subplots(nrows,ncols,figsize=(24,19),squeeze=False)
    fig.subplots_adjust(left=.105,right=.985,bottom=.07,top=.88,hspace=.49,wspace=.39)
    x=np.arange(1,last+1)
    for ri,case in enumerate(contract['cases']):
        obs=case['obs_id'];case_label=('digit '+str(case['digit']))if case.get('digit')is not None else case['family']
        for ci,(metric,(title,unit))in enumerate(METRICS.items()):
            ax=axes[ri,ci];finite_values=[];has_zero_or_negative=False;messages=[]
            ax.axvspan(trained,last+.5,color='#f4e9d4',alpha=.55,zorder=0)
            ax.axvline(trained,color='black',linestyle=':',linewidth=1,zorder=1)
            for mode in mapping:
                table=datasets[mode];values=[];invalid=[];missing=[]
                for iteration in x:
                    raw=table.get((obs,int(iteration)),{}).get(metric);value,status=parsed(raw)
                    values.append(np.nan if value is None else value)
                    if value is not None:
                        finite_values.append(value);has_zero_or_negative|=value<=0
                    elif status=='NOT_RECORDED':missing.append(iteration)
                    else:invalid.append(iteration)
                ax.plot(x,values,color=COLORS[mode],marker='o',markersize=2.7,linewidth=1.2,label=LABELS[mode])
                if invalid:
                    ax.plot(invalid,[.97]*len(invalid),transform=blended_transform_factory(ax.transData,ax.transAxes),linestyle='',marker='x',color=COLORS[mode],markersize=5)
                    messages.append(LABELS[mode]+': nonfinite/nonnumeric '+','.join(map(str,invalid)))
                if missing:
                    ax.plot(missing,[.025]*len(missing),transform=blended_transform_factory(ax.transData,ax.transAxes),linestyle='',marker='|',color=COLORS[mode],markersize=5)
                    messages.append(LABELS[mode]+': missing '+','.join(map(str,missing)))
            positive=[v for v in finite_values if v>0]
            if positive and not has_zero_or_negative:ax.set_yscale('log')
            elif positive:
                # Axis display only: retain exact zero/negative recorded values;
                # this scale is never used as a scientific acceptance threshold.
                ax.set_yscale('symlog',linthresh=min(positive)/10)
            if not finite_values:ax.text(.5,.5,'NOT_RECORDED',transform=ax.transAxes,ha='center',color='.4')
            if messages:ax.text(.01,.01,'\n'.join(messages),transform=ax.transAxes,fontsize=5.7,va='bottom',bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
            ax.set_xlim(.5,last+.5);ax.set_xticks((1,2,4,8,12,16));ax.tick_params(labelsize=8)
            ax.grid(axis='y',alpha=.20);ax.set_ylabel(unit,fontsize=8)
            if ri==0:ax.set_title(title,fontsize=10.5,pad=10)
            if ri==nrows-1:ax.set_xlabel('recurrent iteration',fontsize=9)
            if ci==0:
                noise='clean'if float(case['noise'])==0 else'noise '+str(case['noise'])
                ax.text(-.64,.5,case_label+' / '+noise+'\n'+case['parent_id']+' / A='+str(case['amplitude']),transform=ax.transAxes,ha='center',va='center',fontsize=9.5,fontweight='bold')
    handles=[Line2D([],[],color=COLORS[m],lw=2,label=LABELS[m])for m in mapping]
    handles.append(Line2D([],[],color='black',ls=':',label='trained depth 2'))
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.925),ncol=4,frameon=False,fontsize=11)
    fig.suptitle('VAL public-input recurrent-depth diagnostic — NOT independent TEST\nAll six predeclared clean/noisy cases; iteration 2 trained, iterations 3–16 outside training depth',fontsize=16,y=.98)
    fig.text(.55,.942,'Each row is one fixed case; no best-case filtering. Shading marks extrapolation beyond the trained two iterations.',ha='center',fontsize=11)
    fig.text(.55,.033,'Saved trajectories only. “Actual Maxwell residual” compares the physical model with observed Y, not GT. No fixed-point or GT-convergence guarantee.',ha='center',fontsize=11)
    fig.text(.55,.013,'Log axes where positive; exact zeros use a display-only symmetric-log scale. Missing/nonfinite data stay explicit; endpoints retain availability status.',ha='center',fontsize=9.5)
    fig.savefig(output,dpi=170);plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,default=HERE/'results/recurrent_iterations')
    p.add_argument('--runs',nargs=3,metavar='MODE=DIRECTORY',help='Explicit mapping for all three modes; default uses the declared recurrent_*_seed7_cascade2 run directories')
    a=p.parse_args();mapping=dict(DEFAULT_RUNS)
    if a.runs:
        mapping={}
        for value in a.runs:
            if '='not in value:raise ValueError('--runs entries must be MODE=DIRECTORY')
            mode,directory=value.split('=',1)
            if mode in mapping:raise ValueError('Duplicate mode mapping')
            mapping[mode]=directory
        if set(mapping)!=set(DEFAULT_RUNS):raise ValueError('All three predeclared modes must be retained')
        mapping={mode:mapping[mode]for mode in DEFAULT_RUNS}
    contract=json.loads((HERE/'audit/RECURRENT_ITERATION_DIAGNOSTIC.json').read_text())
    if int(contract['trained_iterations'])!=2 or int(contract['iterations'])!=16:raise ValueError('This plotting plan is explicitly trained-depth2 versus diagnostic-depth16')
    try:paths,data,receipts=input_plan(a.input_root,mapping,contract)
    except FileNotFoundError as error:
        print(str(error),file=sys.stderr);return 2
    flat,document=endpoints(mapping,paths,data,receipts,contract)
    output=HERE/'figures/recurrent/recurrent_depth_val.png';output.parent.mkdir(parents=True,exist_ok=True)
    chart(mapping,data,receipts,contract,output)
    csv_path=a.input_root/'recurrent_depth_endpoints.csv'
    with csv_path.open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    json_path=a.input_root/'recurrent_depth_endpoints.json'
    json_path.write_text(json.dumps(document,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=document['status'],figure=str(output),endpoint_csv=str(csv_path),endpoint_json=str(json_path),model_cases=len(document['endpoints']),endpoint_metric_rows=len(flat),numerical_calls=0),indent=2))
    return 0

if __name__=='__main__':raise SystemExit(main())
