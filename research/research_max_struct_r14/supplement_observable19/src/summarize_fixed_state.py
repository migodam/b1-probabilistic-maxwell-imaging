#!/usr/bin/env python3
"""Summarize archived same-state diagnostics; never loads truth/NPZ or calls physics."""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import csv
import json
from pathlib import Path
import time
import traceback
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]
ARMS=['g19_internal_galerkin','g19_observation_galerkin','g19_anchor_observation_galerkin',
      'g19_internal_right','g19_observation_right','g19_anchor_observation_right']
LABELS=['G basis / PGP (Galerkin)','S basis / PGP (Galerkin)','O0 basis / PGP (Galerkin)',
        'G basis / GP (right-only)','S basis / GP (right-only)','O0 basis / GP (right-only)']
BASIS=['G basis','S basis','O0 basis']
NOTICE='OFFLINE_DIAGNOSTIC — common states from G-PGP trajectory; not recovery predictor'
FIELDS=['case_id','readout_call','arm','status','diagnostic_kind','physical_Phi_before',
        'raw_projected_current_energy_ratio','transported_current_omission_relative_norm',
        'Green_input_omission_relative_norm','PGP_extra_output_trim_relative_to_GPe',
        'full_gradient_to_physical_gradient_angle_degrees',
        'const_gradient_to_physical_gradient_angle_degrees',
        'metric_gradient_to_physical_gradient_angle_degrees',
        'diagnostic_proposal_Phi','diagnostic_Phi_change',
        'diagnostic_step_to_negative_physical_gradient_angle_degrees',
        'diagnostic_material_step_L2','diagnostic_proposal_legal',
        'repeated_internal_readout_saved_m_after_Linf_difference',
        'repeated_internal_readout_saved_m_after_L2_difference','warnings','errors']


def dump(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')


def run(root):
    cpu,wall=time.process_time(),time.perf_counter()
    out=root/'analysis/fixed_state_summary';out.mkdir(parents=True,exist_ok=True)
    receipt=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
                 scope=NOTICE,science_executed=False,truth_arrays_loaded=False,
                 new_SVD_rank_physics_scoring=False,BLAS_threads=1,errors=[])
    try:
        os.environ['MPLCONFIGDIR']=str(out/'.matplotlib')
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from threadpoolctl import threadpool_info
        receipt['numeric_threadpools']=threadpool_info()
        source=root/'audit/fixed_state_v2'
        diagnostics=json.loads((source/'SAME_STATE_SIX_ARM_DIAGNOSTICS.json').read_text())
        common=json.loads((source/'COMMON_STATES.json').read_text())
        states=sorted((r['case_id'],int(r['readout_call'])) for r in common)
        if len(states)!=21 or len(set(states))!=21: raise ValueError('Expected 21 unique common states')
        cases=sorted(set(c for c,k in states))
        if len(cases)!=7 or any(sorted(k for c,k in states if c==case)!=[1,128,512] for case in cases):
            raise ValueError('Expected seven cases x calls1/128/512')
        by={(r['case_id'],int(r['readout_call']),r['arm']):r for r in diagnostics}
        expected={(c,k,a) for c,k in states for a in ARMS}
        if len(diagnostics)!=126 or set(by)!=expected: raise ValueError('Expected 126 unique complete state/arm records')
        rows=[by[c,k,a] for c,k in states for a in ARMS]
        if any(r['status']!='COMPLETE' for r in rows): raise ValueError('Non-COMPLETE archived row')
        common_by={(r['case_id'],int(r['readout_call'])):r for r in common}
        before_delta=max(abs(r['physical_Phi_before']-common_by[r['case_id'],int(r['readout_call'])]['physical_Phi_before']) for r in rows)
        repeat=[by[c,k,ARMS[0]] for c,k in states]
        for r in repeat:
            if r.get('repeated_internal_readout_saved_m_after_Linf_difference') is None: raise ValueError('Missing repeated G-PGP difference')
        repeat_max=max(r['repeated_internal_readout_saved_m_after_Linf_difference'] for r in repeat)
        repeat_L2=max(r['repeated_internal_readout_saved_m_after_L2_difference'] for r in repeat)
        pair_differences={field:max(abs(by[c,k,ARMS[j]][field]-by[c,k,ARMS[j+3]][field]) for c,k in states for j in range(3))
                          for field in ['raw_projected_current_energy_ratio','transported_current_omission_relative_norm']}
        def values(field,arm):
            return [by[c,k,arm].get(field) for c,k in states]
        def stats(values):
            finite=[float(v) for v in values if v is not None and np.isfinite(v)]
            return dict(count=len(finite),missing=len(values)-len(finite),minimum=min(finite) if finite else None,
                        median=float(np.median(finite)) if finite else None,mean=float(np.mean(finite)) if finite else None,
                        maximum=max(finite) if finite else None)
        summary=dict(scope=NOTICE,source_files=[str((source/n).relative_to(root)) for n in ['SAME_STATE_SIX_ARM_DIAGNOSTICS.json','COMMON_STATES.json']],
                     states=21,cases=cases,calls=[1,128,512],arms=[dict(key=k,label=l) for k,l in zip(ARMS,LABELS)],diagnostic_rows=126,
                     common_Phi_before_max_absolute_difference=before_delta,
                     repeated_G_PGP=dict(records=21,reported_m_after_Linf_max_difference=repeat_max,
                                         reported_m_after_L2_max_difference=repeat_L2,
                                         all_reported_differences_exactly_zero=repeat_max==0 and repeat_L2==0,
                                         scope='Checked archived JSON values only; original checkpoint arrays were not reloaded'),
                     paired_basis_projection_metric_max_difference=pair_differences,
                     arm_statistics={a:{f:stats(values(f,a)) for f in ['raw_projected_current_energy_ratio','transported_current_omission_relative_norm','diagnostic_Phi_change','diagnostic_step_to_negative_physical_gradient_angle_degrees']} for a in ARMS},
                     transported_ratio_greater_than_one=[dict(case_id=r['case_id'],readout_call=r['readout_call'],arm=r['arm'],value=r['transported_current_omission_relative_norm']) for r in rows if r.get('transported_current_omission_relative_norm') is not None and r['transported_current_omission_relative_norm']>1],
                     archived_errors=[dict(case_id=r['case_id'],readout_call=r['readout_call'],arm=r['arm'],errors=r['errors']) for r in rows if r.get('errors')],
                     archived_warnings=[dict(case_id=r['case_id'],readout_call=r['readout_call'],arm=r['arm'],warnings=r['warnings']) for r in rows if r.get('warnings')])
        with (out/'DIAGNOSTIC_COMPACT.csv').open('w',newline='',encoding='utf8') as f:
            writer=csv.DictWriter(f,fieldnames=FIELDS);writer.writeheader()
            for r in rows:writer.writerow({k:json.dumps(r.get(k),ensure_ascii=False) if isinstance(r.get(k),(dict,list)) else r.get(k) for k in FIELDS})
        colors=['#2166ac','#b35806','#1b7837','#67a9cf','#e08214','#5aae61']
        figures=[]
        def save(fig,name):
            fig.text(.5,.012,NOTICE,ha='center',va='bottom',fontsize=9)
            for ext in ['png','pdf']:
                p=out/(name+'.'+ext);fig.savefig(p,dpi=160,bbox_inches='tight');figures.append(dict(path=p.name,bytes=p.stat().st_size))
            plt.close(fig)
        # Specific requested object: all three saved readout calls, both projections.
        fig,axes=plt.subplots(3,4,figsize=(21,12))
        for i,call in enumerate([1,128,512]):
            rr=[by['r14_new_0021',call,a] for a in ARMS]
            for j,(field,title,chosen) in enumerate([
                ('raw_projected_current_energy_ratio','Current energy ||Pe||² / ||e||²',range(3)),
                ('transported_current_omission_relative_norm','Transported omission norm ratio',range(3)),
                ('diagnostic_Phi_change','Actual diagnostic Phi(after) − Phi(before)',range(6)),
                ('diagnostic_step_to_negative_physical_gradient_angle_degrees','Step angle to negative physical gradient (degrees)',range(6))]):
                ax=axes[i,j];chosen=list(chosen)
                ax.bar(range(len(chosen)),[rr[q].get(field,np.nan) for q in chosen],color=[colors[q] for q in chosen])
                ax.set_xticks(range(len(chosen)),[BASIS[q] if j<2 else LABELS[q] for q in chosen],rotation=30,ha='right',fontsize=8)
                ax.set_title(f'call {call}: {title}',fontsize=10);ax.grid(axis='y',alpha=.2)
                ax.axhline(0,color='gray',lw=.6)
                if j in [0,1]:ax.axhline(1,color='gray',ls=':',lw=.8)
                if j==3:ax.axhline(90,color='gray',ls=':',lw=.8)
                ax.margins(y=.18)
        fig.suptitle('r14_new_0021 — archived common-input six-arm diagnostics\nEnergy/omission identical across the paired PGP/GP rows; paired steps shown separately',fontsize=14)
        fig.tight_layout(rect=(0,.04,1,.94));save(fig,'R14_NEW_0021_COMMON_STATE_DIAGNOSTICS')
        # Every state and arm appears; signed continuous changes, no recovery labels.
        matrix=np.array([[by[c,k,a].get('diagnostic_Phi_change',np.nan) for a in ARMS] for c,k in states],float)
        bound=float(np.nanmax(abs(matrix))) or 1.
        fig,ax=plt.subplots(figsize=(13,13));im=ax.imshow(matrix,cmap='RdBu_r',vmin=-bound,vmax=bound,aspect='auto')
        ax.set_xticks(range(6),LABELS,rotation=25,ha='right',fontsize=10)
        labels=[f'{c} | call {k}' for c,k in states];ax.set_yticks(range(21),labels,fontsize=9)
        for i in range(21):
            for j in range(6):ax.text(j,i,f'{matrix[i,j]:+.4g}',ha='center',va='center',fontsize=8,color='white' if abs(matrix[i,j])>.55*bound else 'black')
        for y in [2.5,5.5,8.5,11.5,14.5,17.5]:ax.axhline(y,color='black',lw=.5)
        fig.colorbar(im,ax=ax,label='Actual diagnostic Phi(after) − Phi(before)')
        ax.set_title('All 21 common states × six arms: actual archived Phi changes')
        fig.tight_layout(rect=(0,.04,1,1));save(fig,'ALL21_SIX_ARM_PHI_CHANGES')
        # Full-range norm ratios: no truncation to [0,1] and no percent-coverage label.
        fig,axes=plt.subplots(2,1,figsize=(19,11),sharex=True)
        for ax,field,title in zip(axes,['raw_projected_current_energy_ratio','transported_current_omission_relative_norm'],
                                 ['Current projected energy ratio ||Pe||² / ||e||²','Transported omission norm ratio (unclipped; may exceed 1)']):
            for j in range(3):ax.plot(range(21),values(field,ARMS[j]),'o-',color=colors[j],label=BASIS[j]+' (same for PGP/GP)',markersize=4)
            ax.axhline(1,color='gray',ls=':',label='reference ratio = 1');ax.set_title(title);ax.grid(alpha=.2);ax.margins(y=.12);ax.legend(fontsize=9)
        axes[-1].set_xticks(range(21),labels,rotation=65,ha='right',fontsize=8)
        fig.suptitle('All archived common-input states — G, S and O0 basis diagnostics')
        fig.tight_layout(rect=(0,.04,1,.96));save(fig,'ALL21_CURRENT_ENERGY_AND_TRANSPORTED_OMISSION')
        summary['figures']=figures;dump(out/'SUMMARY.json',summary)
        md=['# Archived common-state diagnostic summary','',NOTICE,'',
            'Mechanical summary only: 21 common states = 7 cases × calls 1/128/512; six arms at each state (126 rows). No truth arrays, new physics, scoring, SVD or rank selection.',
            '',f"Reported repeated G-PGP m_after max differences: Linf={repeat_max:.17g}; L2={repeat_L2:.17g} across all 21 states. These values were checked in the archived JSON; original NPZ states were not reloaded.",
            '',f"Common Phi_before maximum difference against COMMON_STATES.json: {before_delta:.17g}.",
            '',f"Transported omission ratios greater than 1: {len(summary['transported_ratio_greater_than_one'])} arm rows (paired projections counted separately). All plotted at their actual values; this norm ratio is not labeled a material-information percentage.",
            '',f"Archived errors: {len(summary['archived_errors'])}; archived warnings: {len(summary['archived_warnings'])}.",'',
            '| Arm | Phi change min | Phi change median | Phi change max | Transported omission max |','|---|---:|---:|---:|---:|']
        for a,l in zip(ARMS,LABELS):
            p=summary['arm_statistics'][a]['diagnostic_Phi_change'];t=summary['arm_statistics'][a]['transported_current_omission_relative_norm']
            md.append(f"| {l} | {p['minimum']:.8g} | {p['median']:.8g} | {p['maximum']:.8g} | {t['maximum']:.8g} |")
        md+=['','All seven case IDs: '+', '.join(cases)+'.','',
             'The object-specific figure contains all three requested calls for r14_new_0021. The other two figures include every common state. Original arm keys and all rows are in DIAGNOSTIC_COMPACT.csv. Scientific interpretation remains with the parent review.','']
        for f in figures:
            if f['path'].endswith('.png'):md.append('!['+f['path']+']('+f['path']+')\n')
        (out/'SUMMARY.md').write_text('\n'.join(md),encoding='utf8')
        receipt.update(status='COMPLETE',source_files=summary['source_files'],states=21,arm_rows=126,figures=figures,
                       repeated_G_PGP=summary['repeated_G_PGP'],all21_states_included=True,transported_ratios_clipped=False,
                       archived_error_rows=len(summary['archived_errors']),archived_warning_rows=len(summary['archived_warnings']))
    except BaseException as exc:
        receipt.update(status='FAILED',errors=[repr(exc)],traceback=traceback.format_exc())
        raise
    finally:
        receipt.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                       completed_utc=datetime.now(timezone.utc).isoformat())
        dump(out/'RECEIPT.json',receipt)
    print(json.dumps({k:receipt[k] for k in ['status','states','arm_rows','CPU_seconds','wall_seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=ROOT)
    run(parser.parse_args().root.resolve())
