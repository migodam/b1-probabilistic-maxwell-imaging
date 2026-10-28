"""Post-seal array-only material, task, spectrum and curvature diagnostics."""
import os
os.environ['MPLBACKEND']='Agg'
from pathlib import Path
import json,csv,time
import numpy as np
from scipy.special import ndtr
import matplotlib.pyplot as plt
from physics import OUT,R4,write
from p3r4.core import MaterialPrior
from material_jets import real_material,decode,material_derivative,tasks


def prior():
    return MaterialPrior(json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text()),
      json.loads((R4/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'])


def main():
    if not (OUT/'E2/PREDICTIONS_SEALED.json').exists():raise RuntimeError('seal first')
    pr=prior();parents=[f'R7_LOCKED_{i:03d}' for i in range(8)]
    records=[];taskrows=[];stability=[];cutoffrows=[];group_rows=[]
    for parent in parents:
        d=OUT/'E2/parents'/parent
        with np.load(d/'BASE_GEOMETRY.npz',allow_pickle=False) as z:g={k:z[k].copy() for k in z.files}
        with np.load(d/'ANCHOR.npz',allow_pickle=False) as z:x=z['x'].copy();m0=z['material'].copy();J=z['J'].copy()
        with np.load(d/'COMPENSATED_JETS.npz',allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
        # Linear coefficient is explicitly NOT cancelled or declared null.
        zz=-np.linalg.lstsq(a['A'],a['c1'],rcond=1e-10)[0]
        U,s,Vh=np.linalg.svd(g['M'],full_matrices=False)
        raw=real_material(a['m2raw'])
        rec=dict(parent=parent,rank=int(g['rank']),material_speed=float(np.linalg.norm(real_material(a['m1']))),
          c1_norm=float(np.linalg.norm(a['c1'])),c2_norm=float(np.linalg.norm(a['c2'])),c3_norm=float(np.linalg.norm(a['c3'])),
          first_order_strong_projectable_fraction=float(np.linalg.norm(a['A']@zz)/np.linalg.norm(a['c1'])),
          decoder2_norm=float(np.linalg.norm(a['decoder2'])),maxwell2_norm=float(np.linalg.norm(a['propagation2'])),
          full2_norm=float(np.linalg.norm(a['b2'])),
          decoder_maxwell2_cosine=float(a['decoder2']@a['propagation2']/max(np.linalg.norm(a['decoder2'])*np.linalg.norm(a['propagation2']),1e-300)),
          material_second_tangent_fraction=float(np.linalg.norm(a['material_coordinate2'])**2/max(np.linalg.norm(raw)**2,1e-300)),
          full_data_min_material_singular=float(g['generalized_joint_singular'][-1]),
          exact_null_or_cubic_claim=False)
        records.append(rec)
        for cutoff in (1e-6,1e-8,1e-10,1e-12):
            keep=s>max(s[0]*cutoff,1e-10)
            T=Vh[keep].T/s[keep][None,:]
            K=J@T
            _,ss,vv=np.linalg.svd(K,full_matrices=False)
            pure=T@vv[-1]
            f=decode(pr,x+.02*pure)[0]
            cutoffrows.append(dict(parent=parent,cutoff=cutoff,material_rank=int(keep.sum()),
               min_response_singular=float(ss[-1]),max_response_singular=float(ss[0]),
               pure_Maxwell_unit_material_latent_norm=float(np.linalg.norm(pure)),
               linear_02_material_actual_RMS=float(np.linalg.norm(f-m0)/24),
               linear_02_prior_increase=float(.5*(np.linalg.norm(x+.02*pure)**2-np.linalg.norm(x)**2))))
        T=g['T'];K=g['K'];mask=pr.mask;right=mask&(pr.points[:,0]>=0);n=int(mask.sum());nr=int(right.sum())
        L=np.zeros((3,1152));L[0,:576]=24*mask/n;L[1,576:]=24*mask/n
        L[2,:576]=48*right*(m0-pr.host).real/nr;L[2,576:]=48*right*(m0-pr.host).imag/nr
        for gn,sl in [('A',slice(0,256)),('B',slice(256,512)),('JOINT',slice(None))]:
            pinv=np.linalg.pinv(K[sl],rcond=1e-10)
            amplification=np.linalg.norm(L@g['Q']@pinv,axis=1)
            taskrows.append(dict(parent=parent,group=gn,task_local_noise_amplifications=amplification.tolist(),
                material_inverse_amplification=float(np.linalg.norm(pinv,2)),
                meaning='linearized inverse sensitivity only; NOT posterior SD/covariance',threshold_event_derivative_claim=False))
        with np.load(d/'B_FULLJ/RMS_0.0200/ENDPOINT.npz',allow_pickle=False) as z:m1=z['material'].copy()
        rr=json.loads((d/'B_FULLJ/RMS_0.0200/RESULT.json').read_text())
        delta=tasks(pr,m1)-tasks(pr,m0)
        dd=rr['whitened_pair_L2']
        tv=2*ndtr(dd/2)-1
        group_rows.append(dict(parent=parent,material_RMS=rr['material_RMS'],whitened_pair_L2=dd,
           tasks_difference=delta.tolist(),two_fixed_Gaussian_worlds_KL=.5*dd*dd,
           two_fixed_Gaussian_worlds_TV=float(tv),
           scalar_two_point_minimax_lowerbounds=(delta**2*(1-tv)/8).tolist(),
           qualification='standard two-point bound for a future fixed-pair experiment, NOT original-prior Bayes risk or selected-data posterior'))
        prev=None
        for node in sorted((d/'paths/B_FULLJ_CHORD').glob('NODE_*.npz')):
            with np.load(node,allow_pickle=False) as z:
                xx=z['x'].copy();MM=real_material(material_derivative(pr,xx));JJ=z['J'].copy()
                uu,ss,vh=np.linalg.svd(MM,full_matrices=False);kk=ss>max(ss[0]*1e-8,1e-10);TT=vh[kk].T/ss[kk][None,:]
                QQ=MM@TT;KK=JJ@TT;_,ks,kv=np.linalg.svd(KK,full_matrices=False)
                W=QQ@kv[-min(4,len(kv)):].T
                Mframe=z['material_frame'].copy()
                # Projector differences are invariant to SVD vector signs/rotations.
                singular=np.linalg.svd(g['Q'].T@Mframe,compute_uv=False)
                base_projector_fro=float(np.sqrt(max(0.,len(g['Q'].T)+Mframe.shape[1]-2*np.sum(singular**2))))
                weak_angle=None if prev is None else float(np.arccos(np.clip(np.linalg.svd(prev.T@W,compute_uv=False)[-1],-1,1)))
                stability.append(dict(parent=parent,node=node.stem,rank=int(kk.sum()),
                  material_projector_Fro_to_anchor=base_projector_fro,
                  weak_max_principal_angle_to_previous=weak_angle,
                  min_response_singular=float(ks[-1]),relative_J_change_to_anchor=float(np.linalg.norm(JJ-J)/np.linalg.norm(J))))
                prev=W
    dest=OUT/'analysis';dest.mkdir(exist_ok=True)
    write(dest/'MATERIAL_PHYSICS.json',dict(rows=records,rank_cutoff_sensitivity=cutoffrows,
        task_local_sensitivity=taskrows,two_world_tasks=group_rows,stability=stability,new_Maxwell_RHS=0,
        coverage='sealed locked8, post-seal descriptive, no method tuning',
        remaining_first_order_profile_qualification='C fixes material velocity; n2/n3 compensated. A^T c1 need not vanish. Not an exact response-aligned profile.'))
    # Fixed gallery: all parents in frozen index order, never GT-selected.
    fig,axs=plt.subplots(8,4,figsize=(10,19))
    for i,parent in enumerate(parents):
        d=OUT/'E2/parents'/parent
        with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:truth=z['material'].copy()
        with np.load(d/'ANCHOR.npz',allow_pickle=False) as z:anchor=z['material'].copy()
        with np.load(d/'B_FULLJ/RMS_0.0200/ENDPOINT.npz',allow_pickle=False) as z:material=z['material'].copy()
        for j,(im,title) in enumerate([(truth.real,'GT (offline)'),(anchor.real,'Measurement-only anchor'),(material.real,'B: material RMS .02'),((material-anchor).real,'B minus anchor')]):
            ax=axs[i,j]
            if j==3:ax.imshow(im.reshape(24,24),origin='lower',cmap='RdBu_r',vmin=-.08,vmax=.08)
            else:ax.imshow(im.reshape(24,24),origin='lower',cmap='viridis',vmin=.1,vmax=1.9)
            ax.set_xticks([]);ax.set_yticks([])
            if i==0:ax.set_title(title,fontsize=10)
            if j==0:ax.set_ylabel(parent[-3:],rotation=0,labelpad=15)
    fig.suptitle('All eight frozen parents; material witnesses are NOT posterior samples',fontsize=12,y=.995)
    fig.text(.5,.006,'Real material: common range 0.1 to 1.9; difference: common range -0.08 to +0.08',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.016,1,.974));fig.savefig(OUT/'figures/LOCKED_MATERIAL_GALLERY.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(13,4))
    for r in records:ax[0].scatter(r['decoder2_norm'],r['maxwell2_norm']);
    ax[0].set(xscale='log',yscale='log',xlabel='Decoder contribution norm',ylabel='Maxwell propagation contribution norm')
    for parent in parents:
        rr=[s for s in stability if s['parent']==parent]
        ax[1].plot(np.linspace(0,1,len(rr)),[s['material_projector_Fro_to_anchor'] for s in rr],'.-')
        ax[2].plot(np.linspace(0,1,len(rr)),[s['relative_J_change_to_anchor'] for s in rr],'.-')
    ax[1].set(xlabel='Path fraction',ylabel='Material tangent projector change (Fro)')
    ax[2].set(xlabel='Path fraction',ylabel='Relative latent J change')
    fig.tight_layout();fig.savefig(OUT/'figures/CURVATURE_AND_J_STABILITY.png',dpi=180);plt.close(fig)
    print(json.dumps(dict(n=8,new_RHS=0,c1_norm_range=[min(r['c1_norm'] for r in records),max(r['c1_norm'] for r in records)],
        min_material_singular_range=[min(r['full_data_min_material_singular'] for r in records),max(r['full_data_min_material_singular'] for r in records)],
        maximum_linear_strong_fraction=max(r['first_order_strong_projectable_fraction'] for r in records)),indent=2))

if __name__=='__main__':main()
