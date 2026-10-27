"""Plot descriptive saved-path target allocation; no truth or physical calls."""
import os
os.environ['MPLBACKEND']='Agg'
import json
from pathlib import Path
import matplotlib.pyplot as plt
TRACK=Path(__file__).resolve().parent


def main():
    folder=TRACK/'stage_target_mismatch';data=json.loads((folder/'SUMMARY.json').read_text());rows=data['summaries'];x=[24,48,96]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(13,4.3))
    values=[[100*r['property_bound_parent_rate'] for r in rows],
            [r['Probit_initial_clipped_property_x_norm_squared']['mean'] for r in rows],
            [r['Probit_RHS']['mean'] for r in rows]]
    labels=['Native endpoint: any property on box boundary (%)',
            'Mean initial property penalty squared (clip only at init)',
            'Mean paid Probit64 complex RHS / parent']
    for ax,y,label in zip(axes,values,labels):
        ax.plot(x,y,'o-',color='#487f9d');ax.set_xticks(x);ax.set_xlabel('Original Native fresh-state cap');ax.set_ylabel(label);ax.grid(alpha=.17)
        for xx,yy in zip(x,y):ax.annotate(f'{yy:.2f}',(xx,yy),xytext=(0,7),textcoords='offset points',ha='center',fontsize=9)
        ax.margins(y=.18)
    fig.suptitle('CDEV64 saved paths: Native target and Probit warm-start allocation')
    fig.text(.015,.015,'All64 parents, same Probit64 cap. Descriptive means/rates; no causal or GT-risk claim; no new solves. Property clip [1e-4,1-1e-4] is initialization only.',fontsize=8)
    fig.tight_layout(rect=(0,.055,1,.93));fig.savefig(folder/'STAGE_TARGET_ALLOCATION.png',dpi=180);fig.savefig(folder/'STAGE_TARGET_ALLOCATION.pdf');plt.close(fig)
    (folder/'FIGURE_SOURCE.json').write_text(json.dumps(dict(source='SUMMARY.json',parent_source='PARENT_STAGE_TARGETS.csv',parents=64,
        figures=['STAGE_TARGET_ALLOCATION.png','STAGE_TARGET_ALLOCATION.pdf'],GT_read=False,new_RHS=0,meaning='descriptive target mismatch/warm start allocation; not proof of causality'),indent=2)+'\n')


if __name__=='__main__':main()
