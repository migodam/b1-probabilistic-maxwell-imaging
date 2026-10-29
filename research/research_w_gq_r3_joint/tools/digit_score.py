"""Post-seal digit8 scoring/figures. No label/template classification or new forward solve."""
from pathlib import Path
import os,json,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy.ndimage import label
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from digit_stress import ROOT,EX,TX,TY,SIGMA,dump

def background_labels(fg):
 lab,n=label(~fg,structure=np.array([[0,1,0],[1,1,1],[0,1,0]]));edge=set(np.r_[lab[0],lab[-1],lab[:,0],lab[:,-1]])-{0};inside=[k for k in range(1,n+1) if k not in edge];return lab,edge,inside

def main():
 out=EX/'results/digit8'
 if not (out/'ONLINE_SEALED.json').exists():raise RuntimeError('No private truth read before seal')
 if not (EX/'DISPLAY_SCORING_RULE.json').exists():raise RuntimeError('Freeze display before scoring')
 s=json.loads((out/'ONLINE_SEALED.json').read_text());a=np.load(out/'endpoint.npz');t=np.load(EX/'private/truth.npz');z=a['z'];zt=t['z'];chi=z[12:61].reshape(7,7);loss=z[61:].reshape(7,7);ct=zt[12:61].reshape(7,7);lt=zt[61:].reshape(7,7);mask=t['mask'];y=np.load(EX/'online/observation.npz')['y'];cx,cy=np.meshgrid((TX[:-1]+TX[1:])/2,(TY[:-1]+TY[1:])/2)
 def center(q):return np.array([(q*cx).sum(),(q*cy).sum()])/max(q.sum(),1e-15)
 lab,edge,ids=background_labels(mask);reps=[]
 for k in ids:
  points=np.argwhere(lab==k);cen=points.mean(axis=0);reps.append(points[np.argmin(np.sum((points-cen)**2,axis=1))])
 checks=[]
 for thresh in [.5,1.,2.]:
  support=chi>=thresh;ll,ee,ii=background_labels(support);which=[int(ll[tuple(p)]) for p in reps];valid=len(which)==2 and all(v>0 and v not in ee for v in which) and which[0]!=which[1]
  checks.append(dict(threshold=thresh,support_IoU=float(np.sum(support&mask)/max(np.sum(support|mask),1)),estimated_enclosed_holes=len(ii),true_hole_representative_components=which,two_true_holes_remain_distinct_enclosed_background=bool(valid)))
 denom=max(float(chi.sum()),1e-15);padding=np.zeros((7,7),bool);padding[:,[0,6]]=True;holes=np.isin(lab,ids);res=(a['data']-y)/SIGMA;cleanres=(a['data']-t['clean'])/SIGMA;wallE=float(np.sqrt(np.mean((z[:6]-zt[:6])**2)));wallL=float(np.sqrt(np.mean((z[6:12]-zt[6:12])**2)))
 metrics=dict(target_complex_contrast_relative_L2=float(np.linalg.norm((chi-ct)+1j*(loss-lt))/np.linalg.norm(ct+1j*lt)),target_epsilon_complex_relative_L2_f1=float(np.linalg.norm((chi-ct)+1j*(loss-lt))/np.linalg.norm(1+ct+1j*lt)),contrast_RMSE=float(np.sqrt(np.mean((chi-ct)**2))),loss_RMSE=float(np.sqrt(np.mean((loss-lt)**2))),contrast_centroid_error=float(np.linalg.norm(center(chi)-center(ct))),offstroke_contrast_mass_fraction=float(chi[~mask].sum()/denom),padding_contrast_mass_fraction=float(chi[padding].sum()/denom),truehole_contrast_mass_fraction=float(chi[holes].sum()/denom),wall_epsilon_RMSE=wallE,wall_loss_RMSE=wallL,observation_white_norm=float(np.linalg.norm(res)),observation_white_RMS=float(np.linalg.norm(res)/np.sqrt(res.size*2)),same_discrete_clean_prediction_white_RMS=float(np.linalg.norm(cleanres)/np.sqrt(cleanres.size*2)),threshold_diagnostics=checks)
 fig,ax=plt.subplots(2,3,figsize=(11,7),layout='constrained')
 for row,arrs,bounds,titles in [(0,[ct,chi,chi-ct],[(0,8),(0,8),(-8,8)],['True target contrast','Reconstructed FREE49-cell contrast','Estimate minus truth']),(1,[lt,loss,chi>=1],[(0,.25),(0,.25),(0,1)],['True independent loss','Reconstructed FREE49-cell loss','Fixed example support: contrast ≥ 1'])]:
  for j,(arr,bound,title) in enumerate(zip(arrs,bounds,titles)):
   im=ax[row,j].imshow(arr,origin='lower',extent=[TX[0],TX[-1],TY[0],TY[-1]],interpolation='nearest',aspect='equal',vmin=bound[0],vmax=bound[1],cmap='coolwarm' if row==0 and j==2 else 'viridis');ax[row,j].set(title=title,xlabel='x / reference wavelength',ylabel='y / reference wavelength');fig.colorbar(im,ax=ax[row,j],shrink=.8)
 fig.suptitle('Digit8 NEW nonparametric stress: no ellipse, template classifier or supplied stroke support\n7×7 free complex target + unknown wall12; cell/stroke0.2λ₀, glyph1.0×1.4λ₀\nSame-discrete exploratory inversion; generated glyphs1..9 are NOT nine reconstructions',fontsize=12);fig.savefig(EX/'truth_vs_reconstruction.png',dpi=190,bbox_inches='tight');fig.savefig(EX/'truth_vs_reconstruction.pdf',bbox_inches='tight');plt.close(fig)
 records=[json.loads(l) for l in (EX/'PHYSICS_LEDGER.jsonl').read_text().splitlines()];cost=dict(full_fe=sum(r['factorizations'] for r in records),rhs=sum(r['rhs'] for r in records));summary=dict(status='ACTUAL_NONPARAMETRIC_EXPLORATORY_RECONSTRUCTION_NOT_PRIMARY_BLIND_PASS',generated_glyph_ids=list(range(1,10)),PDE_generated_digit_ids=[8],reconstructed_digit_ids=[8],not_reconstructed_digit_ids=[1,2,3,4,5,6,7,9],grid_shape_yx=[7,7],target_unknown_complex_cells=49,total_unknown_real_parameters=110,physics='full2D scalarTM P3, original A2 frequencies4TX16RX, sigma.01 per real/imag',same_discrete_truth_inverse=True,independent_finer_validation=False,shape_or_digit_template_online=False,spatial_regularization=False,stops=s['starts'],selected_by=s['selection'],objective=s['objective'],metrics=metrics,cost=cost,primary_R3_gates_unchanged=True,no_NEW_PASS_gate=True,cell_width_lambda0=.2,stroke_width_lambda0=.2,glyph_total_size_lambda0=[1.,1.4],public_target_ROI_lambda0=[-.7,.7,1.4,2.8],notes=['Free target contrast and loss are independently bounded cellwise; no digit/support prior.','Iteration/call/time limits and lack of spatial regularization prevent a failed image being interpreted as physical impossibility.','Two-hole checks at fixed thresholds are descriptive, not selected or calibrated.'])
 dump(EX/'SUMMARY.json',summary);dump(out/'SCORES.json',summary)
 report='# 数字8：新非参数全波探索\n\n只对数字8生成了真实PDE观测并进行了实际110实变量反演；1..9图片是公共字形定义，不能称为九个成像结果。该探索与原椭圆盲测分开。\n\n未知量为墙6个实部/6个loss，以及7×7目标全部49个contrast和49个loss。无模板分类、无ellipse、未提供笔画support；目标为独立box且无空间正则。笔画宽0.2λ₀，对象1.0×1.4λ₀，两个真孔约0.6×0.4λ₀。模型材料格0.2λ₀；FEM公共targetstep0.1λ₀、P3，生成与反演同离散，不是独立细网格PASS。\n\n完整110变量weighted-adjoint L-BFGS-B，每次完整梯度3FE/24RHS，至多两start×40iter，真实calls/450FE/3600RHS/18:30先到停止。初始化为全target统一contrast0.1/loss0.03，不知道数字ID或真墙。原频率、四发十六收和每实虚坐标σ0.01不变。\n\n结果图truth_vs_reconstruction.png的contrast色轴固定0..8、loss固定0..0.25；不得自动拉色轴制造成功。support示例χ≥1；以下同时报告0.5/1/2阈值，不选择最好者，也不新增PASS门。\n\n```json\n'+json.dumps(dict(stops=s['starts'],metrics=metrics,cost=cost),indent=2)+'\n```\n\n以上属于限定时间、粗材料格、单个预声明digit8的实际探索。未收敛/失败不得解释为物理无解；不能证明1..9复杂目标已恢复，也不改变原R3结论。\n'
 (EX/'REPORT_ZH.md').write_text(report);print(json.dumps(dict(status=summary['status'],metrics=metrics,cost=cost)))
if __name__=='__main__':main()
