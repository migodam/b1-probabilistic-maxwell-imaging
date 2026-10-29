from core import *
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import label,binary_fill_holes
from scipy.ndimage import uniform_filter

def structural_similarity(x,y,data_range=1):
 """7x7 uniform-window SSIM, samplecovariance, valid3pixelcroppedinterior."""
 ux=uniform_filter(x,size=7);uy=uniform_filter(y,size=7);covnorm=49/48
 vx=covnorm*(uniform_filter(x*x,size=7)-ux*ux);vy=covnorm*(uniform_filter(y*y,size=7)-uy*uy);vxy=covnorm*(uniform_filter(x*y,size=7)-ux*uy)
 c1=(.01*data_range)**2;c2=(.03*data_range)**2
 sm=((2*ux*uy+c1)*(2*vxy+c2))/((ux*ux+uy*uy+c1)*(vx+vy+c2))
 return float(sm[3:-3,3:-3].mean())
METHODS=['FULL256','AGENT_A','AGENT_B','AGENT_C'];NAMES=['Full256','A exact support','B total field','C global bound']

def topology(a):
 mask=a.reshape(24,24)>=.5
 _,cc=label(mask,structure=np.array([[0,1,0],[1,1,1],[0,1,0]]))
 inv,n=label(~mask,structure=np.ones((3,3)));edge=set(np.r_[inv[0],inv[-1],inv[:,0],inv[:,-1]])
 return int(cc),sum(k not in edge for k in range(1,n+1))

def costs(path):
 rows=[json.loads(z) for z in path.read_text().splitlines() if z.strip()]
 rhs=sum(r.get('RHS',0) for r in rows);lu=sum(r.get('LU',0) for r in rows)
 # Component times inside EVALUATE_COMPLETE are nested: do not sum twice.
 top=[r for r in rows if r.get('event')=='EVALUATE_COMPLETE' or (r.get('method') in ('FULL256','EXACT_SUPPORT') and r.get('wall_seconds') is not None)]
 kernel=sum(r.get('wall_seconds',0) for r in top)
 ev=[r for r in rows if r.get('event')=='EVALUATE_COMPLETE'];forward=[r for r in ev if not r.get('need_jac')]
 return dict(RHS=rhs,LU=lu,corrected_kernel_wall=kernel,trial_evaluations=len(forward),trial_fallbacks=sum(r.get('fallback',False) for r in forward),total_field_evaluations=len(ev),green_model_setup_component=sum(r.get('wall_seconds',0) for r in rows if r.get('event')=='CURRENT_SHAPE_PROJECTED'))

def main():
 if not (P/'audit/BLIND_READ_N.json').exists() or not (P/'audit/BLIND_READ_P.json').exists():raise RuntimeError('blindvisualrecordsrequired')
 key=json.loads((P/'audit/BLIND_KEY_SEALED.json').read_text())['order'];visual={}
 representation=json.loads((P/'audit/DECODER_VISUAL_QUALIFICATION.json').read_text()) if (P/'audit/DECODER_VISUAL_QUALIFICATION.json').exists() else {}
 for ds in ['N','P']:
  read=json.loads((P/'audit'/f'BLIND_READ_{ds}.json').read_text())['panels']
  for j,oid in enumerate(key[ds]):visual[oid]=dict(panel=chr(65+j),**read[chr(65+j)])
 all_rows=[];arrays={};summaries={};ro=load_setup()['ROI'];pr=prior();nlabels=np.arange(1,10)
 for ds in ['N','P']:
  gts=[];rec={m:[] for m in METHODS};mats={m:[] for m in METHODS};gtmat=[];fs=[];ys=[]
  for i in range(9):
   oid=f'{ds}_{i:03d}';gt=dict(np.load(P/'data/truth'/f'{oid}.npz'));gts.append(gt['occupancy']);gtmat.append(gt['material']);ys.append(np.load(P/'data/online'/f'{oid}.npz')['y']);ccgt,hgt=topology(gt['occupancy'])
   for m in METHODS:
    base=P/'results'/ds/m/oid;z=dict(np.load(base/'PREDICTION.npz'));end=json.loads((base/'PREDICTION_SEALED.json').read_text());cost=costs(base/'COST_LEDGER.jsonl');aa=z['occupancy'];chi=z['material'];rec[m].append(aa);mats[m].append(chi);mask=aa>=.5;truthmask=gt['occupancy']>=.5;inter=np.sum(mask&truthmask);union=np.sum(mask|truthmask);cc,h=topology(aa)
    err=chi-gt['material'];occ=aa-gt['occupancy'];val=dict(dataset=ds,digit=i+1,parent=oid,method=m,complex_MSE=float(np.mean(abs(err)**2)),complex_RMS=float(np.sqrt(np.mean(abs(err)**2))),complex_MAE=float(np.mean(abs(err))),ROI_complex_MSE=float(np.mean(abs(err[ro])**2)),real_MSE=float(np.mean(err.real**2)),imag_MSE=float(np.mean(err.imag**2)),occupancy_MSE=float(np.mean(occ**2)),ROI_occupancy_MSE=float(np.mean(occ[ro]**2)),IoU=float(inter/max(union,1)),Dice=float(2*inter/max(mask.sum()+truthmask.sum(),1)),SSIM_full_grid=float(structural_similarity(gt['occupancy'].reshape(24,24),aa.reshape(24,24),data_range=1)),GT_components=ccgt,reconstructed_components=cc,GT_holes=hgt,reconstructed_holes=h,GT_visual_representability_status=('PHYSICS_NATIVE_CLEAR_TARGET' if ds=='N' else representation.get('digits',{}).get(oid,dict(visual='LEGIBLE_THIS_FIT'))['visual']),GT_noise_Phi=phi(gt['Fmean'],ys[-1]),Phi=end['Phi'],white_residual=float(np.sqrt(2*end['Phi'])),noise_ball=end['noise_ball_pass'],status=end['status'],J_calls=end['J_calls'],forward_calls=end['forward_calls'],wall_seconds=end['wall_seconds'],setup_seconds=end['setup_seconds'],raw_kernel_wall=end['physics_kernel_wall'],**cost,final_exact_check_white_error=end['actual_final_white_approx_error'],strict_support=int(np.count_nonzero(aa[ro]>0)),GT_strict_support=int(np.count_nonzero(gt['occupancy'][ro]>0)),delta_real=float(z['delta'].real),delta_imag=float(z['delta'].imag),blind_panel=visual[oid]['panel'],blind_read=visual[oid]['read'],visual_score_0_to_5=visual[oid]['score'],visual_reason=visual[oid]['reason'])
    if ds=='P':val['estimated_latent_norm']=float(np.linalg.norm(z['x'][:64]));val['GT_latent_norm']=float(np.linalg.norm(gt['xgt'][:64]))
    all_rows.append(val)
  gts=np.stack(gts);gm=np.stack(gtmat);ys=np.stack(ys);arrays[ds+'_GT_occupancy']=gts;arrays[ds+'_GT_material']=gm;arrays[ds+'_y']=ys
  for m in METHODS:rec[m]=np.stack(rec[m]);mats[m]=np.stack(mats[m]);arrays[ds+'_'+m+'_occupancy']=rec[m];arrays[ds+'_'+m+'_material']=mats[m];arrays[ds+'_'+m+'_mask']=rec[m]>=.5
  prefix='' if ds=='N' else 'P_'
  for labelname,data in [('GT',gts),('BASELINE_RECON',rec['FULL256']),('AGENT_A_RECON',rec['AGENT_A']),('AGENT_B_RECON',rec['AGENT_B']),('AGENT_C_RECON',rec['AGENT_C'])]:
   fig,axs=plt.subplots(3,3,figsize=(8,9));fig.suptitle(f'{ds}: {labelname}; physical 24x24 grid')
   for i,ax in enumerate(axs.ravel()):ax.imshow(data[i].reshape(24,24),vmin=0,vmax=1,cmap='gray_r',interpolation='nearest');ax.set_title(str(i+1));ax.axis('off')
   fig.tight_layout();fig.savefig(P/'figures'/f'{prefix}DIGITS_{labelname}_1_TO_9.png',dpi=140);plt.close(fig)
  fig,axs=plt.subplots(9,5,figsize=(11,18));datas=[gts]+[rec[m] for m in METHODS];names=['GT']+NAMES
  for i in range(9):
   for j in range(5):
    ax=axs[i,j];ax.imshow(datas[j][i].reshape(24,24),vmin=0,vmax=1,cmap='gray_r',interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
    if i==0:ax.set_title(names[j])
    if j==0:ax.set_ylabel(str(i+1),rotation=0,labelpad=10,fontsize=15)
  fig.suptitle(f'{ds}: same inverse objective, initialization and iteration cap; no templates');fig.tight_layout();fig.savefig(P/'figures'/f'{prefix}DIGITS_SIDE_BY_SIDE.png',dpi=140);plt.close(fig)
  fig,axs=plt.subplots(9,4,figsize=(9,18));names=['GT Re(chi)','Full Re(chi)','GT Im(chi)','Full Im(chi)'];dd=[gm.real,mats['FULL256'].real,gm.imag,mats['FULL256'].imag]
  for i in range(9):
   for j in range(4):
    ax=axs[i,j];ax.imshow(dd[j][i].reshape(24,24),vmin=.1 if j<2 else .04,vmax=1.9 if j<2 else .22,cmap='viridis',interpolation='nearest');ax.axis('off')
    if i==0:ax.set_title(names[j])
  fig.suptitle(f'{ds}: actual complex material');fig.tight_layout();fig.savefig(P/'figures'/f'{prefix}DIGITS_MATERIAL_REAL_IMAG.png',dpi=140);plt.close(fig)
  fig,axs=plt.subplots(3,3,figsize=(9,9));vmax=max(abs(mats['FULL256']-gm).max(),1e-12)
  for i,ax in enumerate(axs.ravel()):im=ax.imshow(abs(mats['FULL256'][i]-gm[i]).reshape(24,24),vmin=0,vmax=vmax,cmap='magma',interpolation='nearest');ax.set_title(str(i+1));ax.axis('off')
  fig.suptitle(f'{ds}: |reconstructed chi - true chi|, common scale');fig.tight_layout();fig.savefig(P/'figures'/f'{prefix}DIGITS_RECONSTRUCTION_ERRORS.png',dpi=140);plt.close(fig)
  rr=[r for r in all_rows if r['dataset']==ds];summary={}
  for m in METHODS:
   r=[q for q in rr if q['method']==m];summary[m]={k:float(np.mean([q[k] for q in r])) for k in ['complex_MSE','ROI_complex_MSE','IoU','Dice','SSIM_full_grid','wall_seconds','corrected_kernel_wall','strict_support']};summary[m].update(noise_ball_pass=sum(q['noise_ball'] for q in r),RHS=sum(q['RHS'] for q in r),LU=sum(q['LU'] for q in r),total_wall_seconds=sum(q['wall_seconds'] for q in r),trial_fallbacks=sum(q['trial_fallbacks'] for q in r),trial_evaluations=sum(q['trial_evaluations'] for q in r),max_material_difference_vs_Full=float(np.max(abs(mats[m]-mats['FULL256']))))
  summaries[ds]=summary
 fig,axs=plt.subplots(1,2,figsize=(12,4.5))
 for ds,ax in zip(['N','P'],axs):
  for m,name in zip(METHODS,NAMES):
   r=[q for q in all_rows if q['dataset']==ds and q['method']==m];ax.scatter([q['wall_seconds'] for q in r],[q['complex_MSE'] for q in r],label=name,alpha=.8)
  ax.set_title(ds+' per-digit end-to-end controlled-job window');ax.set_xlabel('seconds');ax.set_ylabel('full-grid complex MSE');ax.legend(fontsize=8)
 fig.tight_layout();fig.savefig(P/'figures/DIGITS_COST_VS_QUALITY.png',dpi=150);plt.close(fig)
 keys=list(dict.fromkeys(k for r in all_rows for k in r))
 with (P/'DIGIT_PER_NUMBER_RESULTS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(all_rows)
 mrows=[dict(dataset=ds,method=m,**z) for ds,q in summaries.items() for m,z in q.items()]
 keys=list(dict.fromkeys(k for r in mrows for k in r))
 with (P/'METHOD_COMPARISON.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(mrows)
 np.savez_compressed(P/'DIGIT_RECONSTRUCTION_ARRAYS.npz',**arrays)
 write(P/'EXPERIMENT_RESULTS.json',dict(status='NUMERICAL_RESULTS_COMPLETE_RESOURCE_PROTOCOL_FAIL',GT_scoring_utc=time.time(),predictions_sealed_before_scoring=True,visual_read_before_GT_key=True,comparisons=summaries,visual=visual,representation_visual_addendum=representation,statistics_scope='9designeddigitobjects,not populationconfirmation; originalwindowtimingsseparateandpossiblyresourcecontention;not causalspeedtest',all_order_LU_cap_pass=False,SSIM='full-grid localwindows includingknownhost exterior; reportROI MSE separately'))
 print(json.dumps(summaries,indent=2))
if __name__=='__main__':main()
