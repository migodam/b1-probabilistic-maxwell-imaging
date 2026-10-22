"""Final confirmed effects and actual-P1 bracket figures from scored metadata."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write


def main():
    if json.loads((OUT/'MECHANISM_CHECKPOINT.json').read_text())['status']!='COMPLETE':raise RuntimeError('Second cohort must finish')
    dest=OUT/'figures/confirmation';dest.mkdir(parents=True,exist_ok=False)
    result=json.loads((OUT/'L/mechanism_confirmation/evaluation/RESULT.json').read_text())
    names=list(result['family_mean_risk']);values=list(result['family_mean_risk'].values())
    labels=['Full geometry','Anchor only','Capacity/zero','Permuted geometry','Sensitivity only','Backprojection only','Random decoder envelope']
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),constrained_layout=True)
    axes[0].barh(np.arange(7),values);axes[0].set(yticks=np.arange(7),yticklabels=labels,xlabel='Full-grid material MSE',title='Independent second256; mean losses of three seeds');axes[0].invert_yaxis()
    pairs=result['primary_comparisons']
    effects=100*np.array([r['relative_risk_reduction'] for r in pairs]);cis=100*np.array([r['relative_reduction_CI'] for r in pairs])
    axes[1].errorbar(effects,np.arange(3),xerr=np.stack((effects-cis[:,0],cis[:,1]-effects)),fmt='o',capsize=5)
    axes[1].axvline(0,color='black',ls='--');axes[1].set(yticks=np.arange(3),yticklabels=['Full vs backprojection','Full vs sensitivity','Sensitivity vs random envelope'],xlabel='Relative material-risk reduction (%)',title='Three preregistered contrasts; 98.333% paired CIs');axes[1].invert_yaxis()
    for ext in ('png','pdf'):fig.savefig(dest/('MECHANISM_CONFIRMATION.'+ext),dpi=170)
    plt.close(fig)
    bracket=json.loads((OUT/'C/jobs/C_projector_brackets/RESULTS.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10.5,4),constrained_layout=True)
    for parent in bracket['parents']:
        rows=parent['rows'];h=[r['h'] for r in rows]
        axes[0].plot(h,[r['normal_leakage_fraction'] for r in rows],'o-',label=parent['parent'])
        axes[1].plot(h,[r['loop_material_bracket_relative_error'] for r in rows],'o-',label=parent['parent'])
    axes[0].set(xlabel='Finite-difference scale h (material RMS)',ylabel='Normal / total bracket norm',ylim=(0,1),title='Actual P1 rank4 weak distribution')
    axes[1].set(xlabel='h',ylabel='Reverse-loop vs bracket relative error',title='Two fixed scales; no exact-fiber claim')
    for ax in axes:ax.legend()
    for ext in ('png','pdf'):fig.savefig(dest/('P1_FROBENIUS_DIAGNOSTIC.'+ext),dpi=170)
    plt.close(fig)
    write(dest/'MANIFEST.json',dict(status='COMPLETE',new_RHS=0,GT_read=False,
        figures=['MECHANISM_CONFIRMATION.png','P1_FROBENIUS_DIAGNOSTIC.png'],
        risk_scope='256 independent confirmation parents; three fixed-model losses per parent',
        bracket_scope='Two previously exposed actualP1 anchors, two fixed scales, new1920RHS recorded in C; not toy algebra'))


if __name__=='__main__':main()
