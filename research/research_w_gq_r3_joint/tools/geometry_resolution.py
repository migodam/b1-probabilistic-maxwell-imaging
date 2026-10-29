"""Zero-PDE continuous-geometry/raster checks of three already sealed FullGN endpoints."""
from pathlib import Path
import os,json,csv
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/geometry_resolution';OUT.mkdir(exist_ok=True);FIG=ROOT/'figures'
CASES=['blind_normal','blind_cavity','blind_strong'];LABELS=['Normal wall / circle','Boundary-connected air defect\n(historical ID: blind_cavity)','Strong target / ellipse']
def mask(z,x,y):
    dx=x[None,:]-z[12];dy=y[:,None]-z[13];c=np.cos(z[18]);s=np.sin(z[18]);u=c*dx+s*dy;v=-s*dx+c*dy;return (u/z[14])**2+(v/z[15])**2<=1

def contour(z,n=8192):
    # Uniform arclength samples avoid biasing the mean boundary distance toward
    # the slow polar-angle part of an eccentric ellipse.
    a=np.linspace(0,2*np.pi,32769);c=np.cos(z[18]);s=np.sin(z[18]);x=z[12]+c*z[14]*np.cos(a)-s*z[15]*np.sin(a);y=z[13]+s*z[14]*np.cos(a)+c*z[15]*np.sin(a);arc=np.r_[0,np.cumsum(np.hypot(np.diff(x),np.diff(y)))];q=np.arange(n)*arc[-1]/n
    return np.column_stack((np.interp(q,arc,x),np.interp(q,arc,y))),float(arc[-1]/n)

def main():
    rows=[];fig,axes=plt.subplots(1,3,figsize=(14,5.4),layout='constrained')
    for i,case in enumerate(CASES):
        d=ROOT/'results/online'/f'{case}_r0'/'FULL_GN'
        if not (d/'ONLINE_SEALED.json').exists() or not (d/'SCORES.json').exists():raise RuntimeError('Sealed and already scored endpoint required')
        z=np.load(d/'endpoint.npz')['z'];t=np.load(ROOT/'inputs/private'/case/'truth.npz')['z'];score=json.loads((d/'SCORES.json').read_text())['selected'];ious=[];counts=[]
        for nx,ny in [(401,301),(1601,1201)]:
            x=np.linspace(-1.6,1.6,nx);y=np.linspace(.4,2.8,ny);a=mask(z,x,y);b=mask(t,x,y);inter=int(np.sum(a&b));union=int(np.sum(a|b));ious.append(inter/union);counts.append(dict(nx=nx,ny=ny,dx=float(x[1]-x[0]),dy=float(y[1]-y[0]),intersection=inter,union=union,truth_pixels=int(b.sum()),estimated_pixels=int(a.sum())))
        p,step=contour(z);q,tstep=contour(t);ep=cKDTree(q).query(p)[0];tp=cKDTree(p).query(q)[0];aa=np.sort(z[14:16]);ta=np.sort(t[14:16]);absolute=abs(aa-ta);center=float(np.linalg.norm(z[12:14]-t[12:14]));req=float(np.sqrt(t[14]*t[15]));row=dict(parent=case,noise=0,method='FULL_GN',original_score_IoU=score['shape_IoU'],IoU_401x301=ious[0],IoU_1601x1201=ious[1],IoU_raster_change=ious[1]-ious[0],IoU_coarse_reproduction_error=ious[0]-score['shape_IoU'],minor_semiaxis_abs_error=float(absolute[0]),major_semiaxis_abs_error=float(absolute[1]),minor_semiaxis_relative_error=float(absolute[0]/ta[0]),major_semiaxis_relative_error=float(absolute[1]/ta[1]),center_distance=center,true_equivalent_radius=req,center_distance_over_true_equivalent_radius=center/req,true_minimum_diameter=float(2*ta[0]),true_minimum_diameter_pixels_coarse=float(2*ta[0]/.008),true_minimum_diameter_pixels_fine=float(2*ta[0]/.002),boundary_symmetric_Hausdorff_sampled=float(max(ep.max(),tp.max())),boundary_symmetric_mean_arclength=float((ep.mean()+tp.mean())/2),boundary_symmetric_95percentile=float(np.quantile(np.r_[ep,tp],.95)),boundary_sample_arc_step_max=max(step,tstep),boundary_samples_per_ellipse=8192,units='reference wavelength lambda_0',grid_details=counts,truth_sorted_semiaxes=ta.tolist(),estimated_sorted_semiaxes=aa.tolist(),no_PDE=True,no_new_branch_selection=True,no_gate_change=True)
        rows.append(row);ax=axes[i];ax.fill(q[:,0],q[:,1],color='black',alpha=.07);ax.plot(q[:,0],q[:,1],color='black',lw=1.7,label='Continuous truth ellipse');ax.fill(p[:,0],p[:,1],color='#c83636',alpha=.08);ax.plot(p[:,0],p[:,1],color='#c83636',lw=1.7,ls='--',label='Selected continuous estimate');ax.plot(t[12],t[13],marker='+',color='black',ms=9);ax.plot(z[12],z[13],marker='x',color='#c83636',ms=7);ax.plot([t[12],z[12]],[t[13],z[13]],color='#666666',lw=.8,ls=':');cx,cy=(z[12:14]+t[12:14])/2;ax.set(xlim=(cx-.25,cx+.25),ylim=(cy-.25,cy+.25),aspect='equal',xlabel=r'$x/\lambda_0$',ylabel=r'$y/\lambda_0$',title=LABELS[i]+f"\ncenter {center:.4f} λ₀ = {center/req:.2f} true equivalent radii\nIoU coarse {ious[0]:.4f}; fine {ious[1]:.4f}");ax.grid(alpha=.25);ax.plot([cx-.21,cx-.11],[cy-.21,cy-.21],color='black',lw=2);ax.text(cx-.16,cy-.20,'0.1 λ₀',ha='center',va='bottom',fontsize=9)
        ax.text(.02,-.23,f"|Δminor| {absolute[0]:.4f} λ₀; |Δmajor| {absolute[1]:.4f} λ₀\nmean boundary {row['boundary_symmetric_mean_arclength']:.4f} λ₀; max {row['boundary_symmetric_Hausdorff_sampled']:.4f} λ₀",transform=ax.transAxes,va='top',fontsize=9)
    axes[0].legend(fontsize=8,loc='upper right');fig.suptitle('Continuous target-boundary zoom: identical 0.5 λ₀ physical windows and equal aspect ratio\nRender resolution is not measurement resolution; no refit, no reselection, no new PDE',fontsize=13);fig.savefig(FIG/'R0_FULLGN_CONTOUR_ZOOM.png',dpi=200,bbox_inches='tight');fig.savefig(FIG/'R0_FULLGN_CONTOUR_ZOOM.pdf',bbox_inches='tight');plt.close(fig)
    result=dict(rows=rows,physical_forward_fidelity='Existing FINER_VALIDATION reports cover original normal and every registered r0 FullGN three-gate-positive parent. Those field differences are independent of rasterization checks here.',limits=['Increasing the evaluation grid does not add data or refine FEM and never changes the saved continuous estimate.','Boundary distances use arclength samples, with reported step size; they are geometric estimates, not statistical resolution certificates.','No Rayleigh bound, calibrated posterior/coverage or uniform continuous inverse-resolution claim follows from these three point estimates.','Existing position<=.1 gate can be large compared with target radius; center/radius ratios explicitly expose that scale.'],new_physical_solves=0)
    (OUT/'GEOMETRY_RESOLUTION.json').write_text(json.dumps(result,indent=2)+'\n');fields=[k for k in rows[0] if k not in ['grid_details','truth_sorted_semiaxes','estimated_sorted_semiaxes']]
    with (OUT/'GEOMETRY_RESOLUTION.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
    report='# r0 FullGN：几何显示、FEM 与测量分辨率的区分\n\n本检查新增PDE为0，不更改任何端点、分支选择或冻结门。所有轮廓均来自封存的连续参数；三个面板使用相同0.5λ₀物理视窗和等比例坐标。\n\n| Parent | IoU 401×301 | IoU 1601×1201 | 中心距/真值等效半径 | 短/长半轴绝对差 λ₀ | 平均/最大边界距 λ₀ |\n|---|---:|---:|---:|---:|---:|\n'
    for r in rows:report+=f"| {r['parent']} | {r['IoU_401x301']:.6f} | {r['IoU_1601x1201']:.6f} | {r['center_distance_over_true_equivalent_radius']:.4f} | {r['minor_semiaxis_abs_error']:.5f} / {r['major_semiaxis_abs_error']:.5f} | {r['boundary_symmetric_mean_arclength']:.5f} / {r['boundary_symmetric_Hausdorff_sampled']:.5f} |\n"
    report+='\n**栅格化误差**：原评分步长为0.008λ₀，细评分为0.002λ₀，物理范围相同。二者差异只衡量轮廓采样/显示，不增加测量信息，也不改变FEM。\n\n**FEM误差**：另见各选定端点的FINER_VALIDATION.json，比较同一连续参数在不同公共网格上的完整场预测。该检查不是这里的像素细化，也不构成连续PDE证书。\n\n**测量/反问题分辨率**：要比较中心距与目标自身半径、连续边界和材料误差。位置通过0.1λ₀门不表示目标轮廓已精确对齐。本轮只有三个parent及有限噪声复测，不能从高分辨率画图或单点IoU推出统计分辨率、全局唯一性或校准区间。\n'
    (OUT/'RESOLUTION_NOTE_ZH.md').write_text(report);print(json.dumps([dict(parent=r['parent'],IoU_coarse=r['IoU_401x301'],IoU_fine=r['IoU_1601x1201'],center_radius=r['center_distance_over_true_equivalent_radius'],mean_boundary=r['boundary_symmetric_mean_arclength'],hausdorff=r['boundary_symmetric_Hausdorff_sampled']) for r in rows]))
if __name__=='__main__':main()
