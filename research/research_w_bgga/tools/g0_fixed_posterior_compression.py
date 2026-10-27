"""Offline fixed-posterior response compression, using only saved finite15 data.

This is not a cheap online safety rule, not a likelihood substitution, and not
a target-position Bayes-risk certificate. The finite15 physical table is paid.
"""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'results/g0_fixed_posterior';out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for meta in json.loads((ROOT/'results/g0_stress/RESULTS.json').read_text()):
        case=meta['case']
        with np.load(ROOT/'results/g0_stress'/(case+'.npz'),allow_pickle=False) as z:
            p=z['posterior_weights'];F=z['prior_means'];atoms=z['prior_atoms']
        groups=np.arange(15)//3;pos=np.arange(15)%3
        mass=np.bincount(groups,weights=p,minlength=5);representative=int(np.argmax(mass))
        chosen=F[3*representative+pos]
        selection=np.array([8]) if meta['family']=='C3' else np.arange(32)
        a=F.reshape(15,32)[:,selection];b=chosen.reshape(15,32)[:,selection]
        x=np.concatenate((a.real,a.imag),axis=1)/.01;y=np.concatenate((b.real,b.imag),axis=1)/.01
        distortion=float(p@np.sum((x-y)**2,axis=1))
        mean_shift=float(np.sum((p@(x-y))**2))
        rows.append(dict(case=case,family=meta['family'],physical_parent=meta['physical_parent'],
            wall_variance_fraction=meta['finite_posterior']['wall_variance_fraction_of_prior'],
            selected_wall_index=representative,selected_wall_probability=float(mass[representative]),
            expected_response_distortion=distortion,squared_response_prediction_mean_shift=mean_shift,
            jensen_gap=distortion-mean_shift,wall_point_minimum_joint_TV=1-float(mass.max()),
            task='current conditioned 3-pose response task; not the whole25 fixed-probe family',
            identity='OFFLINE_FIXED_POSTERIOR_PREDICTION_ONLY; original joint weights unchanged',
            not_certified='new likelihood/posterior; target position inference; uniform physical error; cheap online action'))
        np.savez_compressed(out/(case+'.npz'),joint_weights=p,original_response=a,compressed_response=b,
            representative_material=atoms[3*representative,:2],matched_pose_index=pos)
    (out/'RESULTS.json').write_text(json.dumps(rows,indent=2)+'\n')
    with (out/'TABLE.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
    colors={'C1':'#426f9d','C2':'#198b76','C3':'#ad7540','C4':'#ba4d55'}
    for r in rows:
        ax.scatter(r['wall_variance_fraction'],r['expected_response_distortion'],color=colors[r['family']])
        ax.annotate(r['case'],(r['wall_variance_fraction'],r['expected_response_distortion']),fontsize=8,xytext=(3,3),textcoords='offset points')
    ax.set_yscale('symlog',linthresh=1e-5)
    ax.set(xlabel='Finite15 wall variance / prior variance',ylabel='Posterior expected response distortion (whitened squared)',
        title='Fixed-posterior diagnostic only; not likelihood-replacement safety')
    fig.savefig(ROOT/'figures/FIXED_POSTERIOR_TASK_DISTORTION.png',dpi=170)
    fig.savefig(ROOT/'figures/FIXED_POSTERIOR_TASK_DISTORTION.pdf');plt.close(fig)
    print(json.dumps(rows))
if __name__=='__main__':main()
