"""Small delivery figures from sealed historical arrays; no inversion."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time
import numpy as np
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
O=Path(__file__).resolve().parent;STUDY=O.parents[1];REPO=STUDY.parents[1];R9=REPO.parent/'r9_worktree';R8=REPO.parent/'r8_worktree'
sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.energy_diagnostics import Geometry,physical_current,energy
CASES=('obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3')
def read(p):return list(csv.DictReader(p.open()))
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def write(name,rows):
    with (O/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def main():
    start=time.perf_counter();cpu=time.process_time();current=[];newbase=[]
    old=read(STUDY/'delegated/energy_audit/development/single_operator_energy.csv')
    new=read(O/'NEW_CANDIDATE_GREEN_ENERGY.csv');bands=[];floor_data=read(STUDY/'delegated/energy_audit/development/actual_FFT_TSom_affine.csv')
    G=S=None;geom=None;jgt={};pubs={}
    for case in CASES:
        pub=load(R9/'research/research_max_struct_r9/data_public/development/operators.npz');pub.update(load(R9/'research/research_max_struct_r9/data_public/development'/f'{case}.npz'));pubs[case]=pub
        if geom is None:G,S=pub['G'],pub['S'];geom=Geometry(G,S,joint=False)
        gt=load(R9/'research/research_max_struct_r9/sealed_truth/development'/f'{case}.npz')['m']
        jgt[case]=physical_current(G,pub['E0'],gt)
        z=load(R8/'research/research_max_struct_r6/RESULTS/primary'/f'{case}__cie_single'/'prediction.npz')
        for quantity,x,family in [('material',z['m']-gt,'complex'),('true_current',z['j']-jgt[case],'complex')]:
            rr,_,_=geom.decompose(x,family)
            bands.extend(dict(case_id=case,method='R6 CIE ref',quantity=quantity,**r) for r in rr if r['tau']==.1)
        for r in old:
            if r['case_id']==case and r['rule']=='fixed_001' and r['state']=='prediction' and r['tau']=='0.1' and r['quantity'] in ('material','true_current') and r['arm'] in ('cie_dual','global_dual','a_plus'):
                bands.append(dict(r,method={'cie_dual':'R9 Dual','global_dual':'R9 Global','a_plus':'R9 A+'}[r['arm']]))
        for r in new:
            if r['case_id']==case and r['state']=='prediction' and r['tau']=='0.1' and r['quantity'] in ('material','true_current'):
                bands.append(dict(r,method={'R10_twofold':'R10 TF28','R10_green_equal_rank':'R10 G28','R10_green_strong':'R10 G19'}[r['method']]))
        dirs={
            'R9 Global':R9/'research/research_max_struct_r9/RESULTS/development'/f'{case}__global_dual__fixed_001',
            'R9 A+':R9/'research/research_max_struct_r9/RESULTS/development'/f'{case}__a_plus__fixed_001',
            'R10 TF28':STUDY/'RESULTS/denominator_dev'/f'{case}__twofold__galerkin',
            'R10 G19':STUDY/'RESULTS/denominator_dev'/f'{case}__green_strong__galerkin'}
        for name,d in dirs.items():
            for st,mf in enumerate((3,6,10,12)):
                aux=load(d/f'stage_{st:02d}_joint.npz');phys=load(d/f'stage_{st:02d}_physical.npz')
                f=next(r for r in floor_data if r['case_id']==case and r['arm']=='global_dual' and r['rule']=='fixed_001' and r['state']==f'stage_{st:02d}_joint')
                current.append(dict(case_id=case,method=name,index=st,mf=mf,state='joint_stage',auxiliary_current_error_energy=energy(aux['J']-jgt[case]),physical_current_error_energy=energy(phys['j']-jgt[case]),affine_GT_floor_energy=float(f['affine_GT_irreducible_energy'])))
            z=load(d/'prediction.npz');current.append(dict(case_id=case,method=name,index=4,mf=None,state='final_physical_outside_affine_constraint',auxiliary_current_error_energy=None,physical_current_error_energy=energy(z['j']-jgt[case]),affine_GT_floor_energy=None))
    write('FIGURE_BAND_ENERGY_DATA.csv',bands);write('FIGURE_AFFINE_FLOOR_DATA.csv',current)
    methods=['R6 CIE ref','R9 Dual','R9 Global','R9 A+','R10 TF28','R10 G28','R10 G19'];names=('strong','weak','complement');colors=('#c54237','#e9ac3a','#3e79b0')
    for case in CASES:
        fig,axs=plt.subplots(2,4,figsize=(17,7.8),constrained_layout=True)
        for qi,quantity in enumerate(('material','true_current')):
            for oi,op in enumerate(('G_in','G_ob')):
                rr=[r for r in bands if r['case_id']==case and r['quantity']==quantity and r['operator']==op]
                absolute=axs[qi,2*oi];fraction=axs[qi,2*oi+1];ys=np.arange(len(methods))
                positive=[]
                for bi,(band,color) in enumerate(zip(names,colors)):
                    ee=np.array([float(next(r for r in rr if r['method']==m and r['band']==band)['component_energy']) for m in methods]);positive.extend(ee[ee>0].tolist())
                    absolute.barh(ys+(bi-1)*.23,np.where(ee>0,ee,np.nan),height=.22,color=color)
                absolute.set_xscale('log');absolute.set_yticks(ys,methods);absolute.invert_yaxis();absolute.grid(axis='x',alpha=.2);absolute.set_xlabel('absolute squared error energy (log)')
                absolute.set_xlim(max(min(positive)*.3,1e-8),max(positive)*3)
                left=np.zeros(len(methods))
                for band,color in zip(names,colors):
                    pct=np.array([100*float(next(r for r in rr if r['method']==m and r['band']==band)['fraction']) for m in methods])
                    fraction.barh(ys,pct,left=left,height=.7,color=color);left+=pct
                fraction.set_xlim(0,100);fraction.set_yticks(ys,methods);fraction.invert_yaxis();fraction.set_xlabel('fraction of total error energy (%)');fraction.grid(axis='x',alpha=.2)
                title=('material coordinate error' if quantity=='material' else 'physical current error')+' | '+op
                absolute.set_title(title+' | absolute',fontsize=10);fraction.set_title(title+' | percentage',fontsize=10)
                if op=='G_in':absolute.text(.02,.01,'complement rank=0; exact zero omitted on log',transform=absolute.transAxes,fontsize=7)
        fig.legend([Patch(color=c) for c in colors],names,loc='outside lower center',ncol=3)
        fig.suptitle(case[4:8]+' | fixed tau=.1 | final saved physical states | disclosed historical DEV',fontsize=13)
        fig.savefig(O/f'ENERGY_{case[4:8]}_ABSOLUTE_PERCENT.png',dpi=170);fig.savefig(O/f'ENERGY_{case[4:8]}_ABSOLUTE_PERCENT.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(13,5.2),constrained_layout=True);palette={'R9 Global':'#a04f37','R9 A+':'#a174b2','R10 TF28':'#227d83','R10 G19':'#dd9b22'}
    for ax,case in zip(axs,CASES):
        for method,c in palette.items():
            rr=[r for r in current if r['case_id']==case and r['method']==method]
            ax.plot([r['index'] for r in rr[:4]],[r['auxiliary_current_error_energy'] for r in rr[:4]],'-o',color=c,markersize=4)
            ax.plot([r['index'] for r in rr],[r['physical_current_error_energy'] for r in rr],'--s',color=c,markersize=3)
        floor=[r for r in current if r['case_id']==case and r['method']=='R9 Global' and r['index']<4]
        ax.plot(range(4),[r['affine_GT_floor_energy'] for r in floor],'-D',color='black',linewidth=2.2,markersize=4)
        ax.set_yscale('log');ax.set_xticks(range(5),['mf3','mf6','mf10','mf12','Final physical']);ax.set_ylabel('squared current error energy (log)');ax.set_title(case[4:8]);ax.grid(alpha=.2)
        ax.axvspan(3.5,4.4,color='#eeeeee');ax.set_xlim(-.15,4.4)
    handles=[Line2D([0],[0],color=c,label=m) for m,c in palette.items()]+[Line2D([0],[0],color='black',marker='D',label='shared affine GT floor (auxiliary only)'),Line2D([0],[0],color='gray',marker='o',label='auxiliary current error'),Line2D([0],[0],color='gray',linestyle='--',marker='s',label='physical current error')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Auxiliary affine floor is not a lower bound on m-derived / final physical current error',fontsize=12)
    fig.savefig(O/'AFFINE_FLOOR_VS_AUXILIARY_PHYSICAL.png',dpi=180);fig.savefig(O/'AFFINE_FLOOR_VS_AUXILIARY_PHYSICAL.pdf');plt.close(fig)
    (O/'FIGURE_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',scope='b20d/40ac final band bars plus four-method sparse stage current chart; no449state rerun',new_inverse_runs=0,new_offline_physical_LU=2,new_offline_complex_RHS=24,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,zero_handling='Exact zero bands omitted on logarithmic absolute-energy axis; percentages retain zero. No energy normalization used in absolute plots.',GT_status='EXPOSED_OFFLINE_DIAGNOSTIC'),indent=2))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
