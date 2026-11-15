"""Data-only plots of the frozen finite-mode trajectories; zero solver calls."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse
import json
from pathlib import Path
import sys
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();cpu,wall=time.process_time(),time.perf_counter();root=args.root.resolve()
    sys.path.insert(0,str(root/'src'));from run_vector_factorial import dump
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    source=root/'audit/finite_mode_trajectory_v1';rows=json.loads((source/'FINITE_PROGRAM_ROWS.json').read_text())
    with np.load(source/'MODE_AND_IDEAL_PREDICTIONS.npz',allow_pickle=False) as z:
        ideal=z['ideal_material_traces'].copy();rhopowers=z['rho01_powers'].copy()
    folder=source/'FIGURES';folder.mkdir(exist_ok=False);pairs=('00','01','10','11')
    colors={'native_PR':'#2768aa','reset_before_every_update':'#d27325'}
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for index,pair in enumerate(pairs):
        ax=axes.ravel()[index]
        for policy in colors:
            for sign in (1,-1):
                values=[row for row in rows if row['pair']==pair and row['memory_policy']==policy and row['sign']==sign]
                ax.plot([0]+[row['step'] for row in values],[1.]+[row['material_gain'] for row in values],
                    color=colors[policy],ls='-' if sign==1 else '--',label=policy+(' +eps' if sign==1 else ' -eps'))
        ax.plot(range(33),np.linalg.norm(ideal[index],axis=1),'k-.',lw=1.4,label='saved ideal T_pair^k h')
        if pair=='01':ax.plot(range(33),rhopowers,color='#777777',ls=':',lw=2,label='saved ideal01 rho^k')
        ax.axhline(1,color='#aaaaaa',lw=.6);ax.set_title('pair'+pair);ax.set_xlabel('native completed updates k')
        ax.set_ylabel('||m_k - m_root|| / .001');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Declared real leading mode of IDEAL01; fixed ±.001, 32 actual PR/readout updates\nFinite-scale program response only; guard/floating-floor effects remain, no Jacobian or asymptotic claim')
    for extension in ('png','pdf'):fig.savefig(folder/('MATERIAL_MODE_GAIN.'+extension),dpi=170,bbox_inches='tight')
    plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(16,7),constrained_layout=True)
    for column,pair in enumerate(pairs):
        for policy in colors:
            for sign in (1,-1):
                values=[row for row in rows if row['pair']==pair and row['memory_policy']==policy and row['sign']==sign]
                for index,key in enumerate(('current_relative_step','field_relative_step')):
                    axes[index,column].plot([row['step'] for row in values],[row[key] for row in values],
                        color=colors[policy],ls='-' if sign==1 else '--',label=policy+(' +eps' if sign==1 else ' -eps'))
        for index in range(2):
            ax=axes[index,column];ax.set_yscale('log');ax.grid(alpha=.2);ax.set_title('pair'+pair);ax.set_xlabel('completed update k')
            if column==0:ax.set_ylabel(('current J' if index==0 else 'field Q')+' relative step, log')
    axes[0,0].legend(fontsize=7)
    fig.suptitle('Current and total-field steps are separate observables; fixed mode/epsilon and original flags\nRelative step = ||x_k-x_(k-1)|| / ||x_(k-1)||; no new Maxwell/optimizer calls')
    for extension in ('png','pdf'):fig.savefig(folder/('CURRENT_AND_FIELD_RELATIVE_STEPS.'+extension),dpi=170,bbox_inches='tight')
    plt.close(fig)
    dump(source/'PLOT_RECEIPT.json',dict(status='COMPLETE',source='existingFINITE_PROGRAM_ROWS and MODE_AND_IDEAL_PREDICTIONS only',
        figures=['MATERIAL_MODE_GAIN','CURRENT_AND_FIELD_RELATIVE_STEPS'],CPU_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-wall,forward_LU=0,solver_RHS=0,optimizer_calls=0,
        zero_current_steps=sum(row['current_relative_step']==0 for row in rows),
        zero_field_steps=sum(row['field_relative_step']==0 for row in rows),predictions_modified=False))


if __name__=='__main__':main()
