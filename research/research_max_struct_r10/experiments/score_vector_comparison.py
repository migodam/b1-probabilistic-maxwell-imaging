"""Postseal scoring/figures of all 18 frozen vector arms; zero new forwards.

Never import generator designs or resample fine material to choose GT. Cross-grid
material and current comparisons use the declared grid12 effective reference;
the fine/coarse clean-port discrepancy is separately reported.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import argparse,csv,json,time,sys
import numpy as np
from scipy.linalg import svd
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
ARMS=('sourceCIE','GlobalDual','AplusGlobal')


def read(path):return json.loads(Path(path).read_text())
def dump(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def table(path,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader();writer.writerows(rows)
def squared(x):return float(np.vdot(x,x).real)
def relative(x,reference):
    n=squared(reference);return float(np.sqrt(squared(x)/n)) if n else None


def material_metrics(m,t):
    delta=m-t;support=m>.5;truth=t>.5;union=np.count_nonzero(support|truth)
    result=dict(chi_relative_L2=relative(delta,t),
        amplitude_normalized_RMSE=float(np.sqrt(np.mean(delta**2))/(t.max()-t.min())),
        paper_Err=float(np.sqrt(np.mean((delta/(1+t))**2))),
        fixed_IoU=float(np.count_nonzero(support&truth)/union) if union else 1.,
        material_legal=bool(np.isfinite(m).all() and m.min()>=-1e-12 and m.max()<=8+1e-12),
        material_min=float(m.min()),material_max=float(m.max()))
    for name,mask in dict(background=t==0,low=(t>0)&(t<=1),middle=(t>1)&(t<=3),high=(t>3)&(t<=8)).items():
        result[name+'_count']=int(mask.sum())
        result[name+'_RMSE']=float(np.sqrt(np.mean(delta[mask]**2))) if mask.any() else None
        result[name+'_bias']=float(np.mean(delta[mask])) if mask.any() else None
    return result


def receiver_basis(config):
    from research.research_max_struct_r10.src.vector_grid import grid_points,external_geometry
    from research.research_max_struct_r1.src.a9_green import dyadic_green
    points=grid_points(12,.6);positions,_=external_geometry(config['acquisition'])
    S=np.empty((180,5184),complex)
    for i in range(0,60,4):
        g=dyadic_green(positions[i:i+4,None]-points[None],2*np.pi)
        S[3*i:3*(i+4)]=.05**3*g.transpose(0,2,1,3).reshape(12,5184)
    _,singular,Vh=svd(S,full_matrices=False,check_finite=False)
    return Vh.conj().T,singular


def projected_energy(error,V,Gv):
    volume=.05**3;total=volume*squared(error)
    ret=volume*squared(V.conj().T@error) if V is not None else None
    top=volume*squared(Gv.conj().T@error)
    return dict(current_error_energy_h3=total,receiver_retained_error_energy=ret,
        receiver_complement_error_energy=total-ret if ret is not None else None,
        receiver_retained_error_percent=100*ret/total if total and ret is not None else None,
        Green_top32_error_energy=top,Green_uncaptured_error_energy=total-top,
        Green_top32_error_percent=100*top/total if total else None,
        Green_scope='top32 vs orthogonal remainder; not full strong/weak/kernel spectrum')


def figures(destination,cohorts,materials,rows,curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    destination.mkdir(parents=True,exist_ok=False)
    norm=Normalize(0,8);cmap=plt.get_cmap('viridis');coordinates=np.linspace(-.3,.3,13)
    X,Y,Z=np.meshgrid(coordinates,coordinates,coordinates,indexing='ij')
    for cohort in cohorts:
        arrays=[materials[cohort]['GT']]+[materials[cohort].get(arm) for arm in ARMS]
        labels=['GT effective chi','CIE','Global-Dual (ablation)','A+ global']
        fig,axes=plt.subplots(3,4,figsize=(13,9),constrained_layout=True)
        for col,(m,label) in enumerate(zip(arrays,labels)):
            for axis in range(3):
                ax=axes[axis,col]
                if m is None:ax.text(.5,.5,'FAILED',ha='center');ax.set_axis_off();continue
                volume=m.reshape(12,12,12);plane=np.take(volume,6,axis=axis)
                ax.imshow(plane.T,origin='lower',extent=(-.3,.3,-.3,.3),norm=norm,cmap=cmap)
                if axis==0:ax.set_title(label)
                if col==0:ax.set_ylabel(('x','y','z')[axis]+' = +0.025 lambda')
        fig.suptitle(cohort+' | fixed central orthogonal slices | all 1728 materials unknown')
        fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes,shrink=.65,label='contrast chi')
        fig.savefig(destination/(cohort+'_slices.png'),dpi=180);fig.savefig(destination/(cohort+'_slices.pdf'));plt.close(fig)
        fig=plt.figure(figsize=(15,4.5),layout='constrained')
        for col,(m,label) in enumerate(zip(arrays,labels),1):
            ax=fig.add_subplot(1,4,col,projection='3d')
            if m is not None:
                v=m.reshape(12,12,12);colors=cmap(norm(v));colors[...,3]=.85
                ax.voxels(X,Y,Z,v>.5,facecolors=colors,edgecolor=None)
            ax.set(xlim=(-.3,.3),ylim=(-.3,.3),zlim=(-.3,.3),title=label)
            ax.set_box_aspect((1,1,1));ax.view_init(24,38)
        fig.suptitle(cohort+' | displayed support chi > 0.5; colors retain contrast')
        fig.savefig(destination/(cohort+'_volume.png'),dpi=180);fig.savefig(destination/(cohort+'_volume.pdf'));plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    for ax,cohort in zip(axes.flat,cohorts):
        for arm in ARMS:
            data=curves.get((cohort,arm),[])
            if data:ax.semilogy([r['accepted'] for r in data],[r['Phi'] for r in data],label=arm)
        ax.set(title=cohort,xlabel='accepted physical-corrector steps',ylabel='original Phi');ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=8);fig.savefig(destination/'PHYSICAL_CONVERGENCE.png',dpi=180)
    fig.savefig(destination/'PHYSICAL_CONVERGENCE.pdf');plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--authorize-postseal',required=True)
    ap.add_argument('--authorization-note',required=True)
    ap.add_argument('--output',type=Path,default=ROOT/'delegated/vector_score/scored')
    ap.add_argument('--figure-dir',type=Path,default=ROOT/'FIGURES/vector_comparison')
    args=ap.parse_args()
    if args.authorize_postseal!='ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION':raise ValueError('Explicit root authorization required')
    wall,cpu=time.perf_counter(),time.process_time()
    tasks=read(ROOT/'configs/vector_all18.json')['tasks']
    if len(tasks)!=18:raise ValueError('This scorer covers exactly the frozen18 primary tasks')
    sealed=[];failed=[];receipts={}
    for task in tasks:
        path=ROOT/'RESULTS/vector'/task['task_id'];r=read(path/'RECEIPT.json');receipts[task['task_id']]=r
        if r['status']=='SEALED_PUBLIC_PREDICTION_NO_SCORE':
            seal=read(path/'SEALED.json')
            if seal['status']!='PREDICTION_SEALED_BEFORE_ROOT_OFFLINE_SCORING':raise ValueError('Invalid public seal')
            if not (path/'prediction.npz').exists():raise ValueError('Missing sealed prediction')
            sealed.append(dict(task=task,receipt=str(path/'RECEIPT.json'),prediction_bytes=(path/'prediction.npz').stat().st_size))
        elif r['status'].startswith('FAILED'):failed.append(dict(task=task,receipt=str(path/'RECEIPT.json'),reason=r.get('failure')))
        else:raise ValueError('Nonterminal vector task: '+task['task_id'])
    args.output.mkdir(parents=True,exist_ok=False)
    dump(args.output/'ALL_PRIMARY_PUBLIC_SEAL.json',dict(authority=args.authorization_note,sealed=sealed,failed=failed,GT_access_so_far=False))
    config=read(ROOT/'configs/vector_actual_contract.json');Vfull,singular=receiver_basis(config)
    with np.load(ROOT/'data_public/vector/operator_spectral/GIN_TOP32_CAPTURE.npz',allow_pickle=False) as z:Gv=z['right_vectors']
    ranks={}
    for path in (ROOT/'RESULTS/vector').glob('*/TERMINALS.json'):
        for c in read(path).get('cases',[]):
            rank=c['receiver']['rank'];key=c['public_id']
            if key in ranks and ranks[key]!=rank:raise ValueError('Inconsistent paid DPCS ranks')
            ranks[key]=rank
    # First private read occurs only after all18 public terminal checks above.
    truth={};materials={};rows=[];state_rows=[];energy=[];curves={};model_error=[]
    cohorts=list(dict.fromkeys(t['public_id'] for t in tasks))
    for cohort in cohorts:
        case=cohort.split('_')[0]
        if case not in truth:
            with np.load(ROOT/'sealed_truth/vector'/case/'grid12.npz',allow_pickle=False) as z:
                truth[case]={key:z[key].copy() for key in ('material','J','Y_clean')}
        materials[cohort]={'GT':truth[case]['material']}
        if 'cross' in cohort:
            with np.load(ROOT/'sealed_truth/vector'/case/'grid18.npz',allow_pickle=False) as z:yfine=z['Y_clean']
            model_error.append(dict(public_id=cohort,clean_grid12_vs_grid18_relative=relative(truth[case]['Y_clean']-yfine,yfine),
                scope='Independent clean port discretization discrepancy; not added into measurement sigma'))
    for task in tasks:
        ident=task['task_id'];cohort=task['public_id'];arm=task['arm'];case=cohort.split('_')[0]
        path=ROOT/'RESULTS/vector'/ident;r=receipts[ident]
        row=dict(task_id=ident,case_id=case,public_id=cohort,arm=arm,status=r['status'],condition='cross18to12' if 'cross' in cohort else 'matched12')
        if r['status'].startswith('FAILED'):rows.append(row);continue
        with np.load(ROOT/'data_public/vector'/cohort/'observations.npz',allow_pickle=False) as z:Y=z['Y']
        sigma=read(ROOT/'data_public/vector'/cohort/'PUBLIC.json')['sigma'];reference=truth[case]
        with np.load(path/'prediction.npz',allow_pickle=False) as z:
            m=z['material'];prediction=z['prediction'];J=z['J']
            row.update(material_metrics(m,reference['material']),physical_relative_residual=relative(prediction-Y,Y),
                physical_Phi=r['physical_Phi'],whitened_port_squared=squared((prediction-Y)/sigma),
                coarse_reference_current_relative_error=relative(J-reference['J'],reference['J']))
            materials[cohort][arm]=m.copy()
            V=Vfull[:,:ranks[cohort]] if cohort in ranks else None
            energy.append(dict(task_id=ident,state='prediction',reference='grid12_same_material_coarse_model_current',
                **projected_energy(J-reference['J'],V,Gv)))
        row.update(PG_inf=r['PG_inf'],nested_inverse_CPU_seconds=r['CPU_seconds'],nested_inverse_wall_seconds=r['wall_seconds'],
            backend=r['backend'],corrector_requests=r['physical_corrector']['requests'],corrector_accepted=r['physical_corrector']['accepted'],
            corrector_stop=r['physical_corrector']['stop'],physical_LU=r['counts'].get('physical_LU'),
            primal_RHS=r['counts'].get('primal_RHS',r['counts'].get('physical_source_RHS')),
            adjoint_RHS=r['counts'].get('adjoint_RHS',r['counts'].get('physical_adjoint_RHS')),
            auxiliary_iterations=sum(s['iterations'] for s in r['stages']),DPCS_rank=ranks.get(cohort))
        rows.append(row)
        for name in ('physical_initial','stage1_physical','stage2_physical','stage3_physical','physical_accepted40','physical_accepted80','prediction'):
            target=path/(name+'.npz')
            if not target.exists():continue
            with np.load(target,allow_pickle=False) as z:
                sm=z['material'];sj=z['J']
                meta=read(path/(name+'.json'))
                state_rows.append(dict(task_id=ident,state=name,**material_metrics(sm,reference['material']),Phi=meta['Phi']))
                if name.startswith('stage'):
                    jp=path/(name.replace('_physical','_joint')+'.npz')
                    with np.load(jp,allow_pickle=False) as aux:
                        if np.array_equal(aux['material'],sm):
                            ep=sj-reference['J'];gap=aux['J']-sj;ea=aux['J']-reference['J'];h3=.05**3
                            energy.append(dict(task_id=ident,state=name.replace('_physical','_joint'),reference='grid12_coarse_model_current',
                                **projected_energy(ea,V,Gv),physical_error_energy=h3*squared(ep),
                                auxiliary_physical_gap_energy=h3*squared(gap),signed_cross_term=2*h3*float(np.vdot(ep,gap).real)))
                        else:energy.append(dict(task_id=ident,state=name,missing_reason='Joint/physical material mismatch; no extra solve'))
        curves[(cohort,arm)]=[z for z in (json.loads(line) for line in (path/'TRACE.jsonl').read_text().splitlines()) if z.get('kind')=='physical_accepted']
    table(args.output/'FINAL_SCORES.csv',rows);table(args.output/'SAVED_STATES.csv',state_rows)
    table(args.output/'CURRENT_ENERGY.csv',energy);table(args.output/'CROSSGRID_MODEL_ERROR.csv',model_error)
    paired=[]
    for cohort in cohorts:
        by={r['arm']:r for r in rows if r['public_id']==cohort}
        for comparator in ('sourceCIE','GlobalDual'):
            a=by['AplusGlobal'];base=by[comparator]
            available='chi_relative_L2' in a and 'chi_relative_L2' in base
            paired.append(dict(public_id=cohort,reference=comparator,
                delta_chi_relative_L2=a['chi_relative_L2']-base['chi_relative_L2'] if available else None,
                delta_fixed_IoU=a['fixed_IoU']-base['fixed_IoU'] if available else None,
                missing_reason=None if available else 'At least one frozen arm failed'))
    table(args.output/'PAIRED_DIFFERENCES.csv',paired)
    figures(args.figure_dir,cohorts,materials,rows,curves)
    receipt=dict(status='COMPLETE',tasks=len(tasks),sealed=len(sealed),failed=len(failed),independent_materials=4,
        crossgrid_cohorts=2,truth_access='POSTSEAL_OFFLINE_SCORING_ONLY',new_forward_calls=0,new_adjoint_calls=0,
        receiver_SVD=1,receiver_S_assembly=1,Green_SVD=0,Green_top32_reused=True,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        cost_scope='Whole scorer after entry, includes operator assembly/economy receiver SVD/arrays/figures; perarm inverse CPU is nested and not total process lifetime',
        confidence='Four independent objects; no inferential success-rate or universal-dominance claim',
        current_reference='grid12; crossgrid current errors are coarse-model references, not fine-field interpolation')
    dump(args.output/'RECEIPT.json',receipt);print(json.dumps(receipt))


if __name__=='__main__':main()
