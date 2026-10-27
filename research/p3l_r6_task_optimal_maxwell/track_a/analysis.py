"""Offline A analysis: saved arrays and receipts only. No solver imports or calls."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
from pathlib import Path
import json,csv,time,resource
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;OUT=HERE/'analysis';OUT.mkdir(exist_ok=True)
def read(p):return json.loads(Path(p).read_text())
def write(name,x):(OUT/name).write_text(json.dumps(x,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def csvwrite(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def cmean(r):return complex(*r['mean_increment'])
def pack(a):return np.concatenate((a.real.transpose(0,2,1),a.imag.transpose(0,2,1)),axis=2).reshape(len(a),-1)
def error_record(a,b,label,kind):
 delta=abs(cmean(a)-cmean(b));mass=max(abs(np.array(a['source_layer_mass'])-b['source_layer_mass']));phase=abs(a['latent_phase_gt_half']-b['latent_phase_gt_half']);event=abs(a['material_real_increment_gt_002']-b['material_real_increment_gt_002']);sd=b['material_ROI_SD'];q=dict(comparison=label,reference_kind=kind,observation=a['observation'],method=a.get('method',f"GH{a.get('order')}"),mean_complex_abs_error=delta,reference_material_SD=sd,mean_error_SD=delta/max(sd,1e-300),mean_abs_gate_1e_5=delta<=1e-5,mean_SD_gate_0p01=delta<=.01*sd,branch_mass_max_abs_error=float(mass),branch_mass_gate_0p002=mass<=.002,latent_phase_abs_error=phase,latent_phase_gate_0p002=phase<=.002,material_event_abs_error=event,material_event_gate_0p002=event<=.002,SD_relative_difference=abs(a['material_ROI_SD']/sd-1),log_evidence_difference=a['log_evidence']-b['log_evidence']);q['joint_frozen_task_gates_pass']=all(q[k] for k in ['mean_abs_gate_1e_5','mean_SD_gate_0p01','branch_mass_gate_0p002','latent_phase_gate_0p002','material_event_gate_0p002']);return q

def main():
 tick=time.perf_counter();results=read(HERE/'RESULTS.json');rows=results['rows'];local=read(HERE/'LOCAL_REFINEMENT_RESULTS.json')['rows'];comp=read(HERE/'COMPENSATION.json');freeze=read(ROOT/'configs/P0_FREEZE.json')['A'];lookup={(r['method'],r['observation']):r for r in rows};refs={r['observation']:r for r in local if r['order']==21};observ=np.load(HERE/'OBSERVATIONS.npz');obs={str(n):(y,float(s)) for n,y,s in zip(observ['names'],observ['observation'],observ['sigma'])}
 # Different paths and global refinement are reported separately, not silently pooled.
 quals=[]
 for a in local:
  if a['order']==15:quals.append(error_record(a,refs[a['observation']],f"{a['observation']}: GH15->GH21",'independent GN-whitened material quadrature refinement Q1/Q2'))
 pairs=[('S1_scaled_12_12','S1_scaled_32_24'),('S1_scaled_24_16','S1_scaled_32_24'),('S1_unscaled_32_24','S1_scaled_32_24'),('S3_source_10_8','S3_source_16_9'),('S3_quotient_12_6','S3_quotient_20_8')]
 for lo,hi in pairs:
  for oname in obs:quals.append(error_record(lookup[lo,oname],lookup[hi,oname],f'{lo}->{hi}','same conditional submodel global-grid Q1/Q2 refinement; no tail certificate'))
 for method in ['S3_source_10_8','S3_source_16_9','S3_quotient_12_6','S3_quotient_20_8']:
  for oname,ref in refs.items():quals.append(error_record(lookup[method,oname],ref,f'{method}->GN_GH21','global conditional S3 versus local GN-fibre Q1/Q2 reference; global tail unbounded'))
 csvwrite('TASK_QUALIFICATION.csv',quals);write('TASK_QUALIFICATION.json',dict(frozen_thresholds=freeze,records=quals,global_reference_certificate=False))
 flat=[]
 for r in rows:
  flat.append(dict(method=r['method'],observation=r['observation'],sigma=r['sigma'],nodes=r['node_count'],incremental_shared_RHS=r['additional_shared_RHS'],mean_real=r['mean_increment'][0],mean_imag=r['mean_increment'][1],material_SD=r['material_ROI_SD'],host_mass=r['source_layer_mass'][0],transition_mass=r['source_layer_mass'][1],pure_mass=r['source_layer_mass'][2],phase_mass=r['latent_phase_gt_half'],material_event_mass=r['material_real_increment_gt_002'],prior_mass=r['prior_mass_quadrature'],log_evidence=r['log_evidence']))
 csvwrite('GLOBAL_GRID_TASKS.csv',flat)
 # Host controlled-noise trajectories, three fixed scales with same innovation.
 host=[];hostslopes=[]
 for method in ['S1_scaled_32_24','S3_source_16_9','S3_quotient_20_8']:
  rr=sorted([r for r in rows if r['method']==method and r['observation'].startswith('HOST')],key=lambda x:x['sigma']);sx=np.array([r['sigma'] for r in rr]);mu=np.array([abs(cmean(r)) for r in rr]);mass=np.array([1-r['source_layer_mass'][0] for r in rr]);sd=np.array([r['material_ROI_SD'] for r in rr]);hostslopes.append(dict(method=method,mean_magnitude_log_slope=float(np.polyfit(np.log(sx),np.log(mu),1)[0]),nonhost_mass_log_slope=float(np.polyfit(np.log(sx),np.log(mass),1)[0]),material_SD_log_slope=float(np.polyfit(np.log(sx),np.log(sd),1)[0]),scope='three conditional controlled observations, same noise innovation; no population Bayes rate'))
  for r in rr:host.append(dict(method=method,sigma_ratio=r['sigma']/obs['HOST_sigma1.0'][1],sigma=r['sigma'],mean_magnitude=abs(cmean(r)),mean_real=r['mean_increment'][0],mean_imag=r['mean_increment'][1],SD=r['material_ROI_SD'],host_mass=r['source_layer_mass'][0],nonhost_mass=1-r['source_layer_mass'][0]))
 csvwrite('HOST_TRAJECTORIES.csv',host);write('HOST_SLOPES.json',hostslopes)
 # Gauge: fit response order, coefficient ratios, and noise-floor compensation.
 trials=comp['finite_trials'];h=np.array([r['h'] for r in trials]);un=np.array([r['uncompensated_white_norm'] for r in trials]);co=np.array([r['compensated_white_norm'] for r in trials]);coef=comp['operators'][0]['c_norm'];gauge=dict(uncompensated_allpoint_log_slope=float(np.polyfit(np.log(h),np.log(un),1)[0]),uncompensated_last3_log_slope=float(np.polyfit(np.log(h[-3:]),np.log(un[-3:]),1)[0]),cubic_coefficient_norm=coef,uncompensated_norm_over_h_cubed=(un/h**3).tolist(),ratio_to_cubic_coefficient=(un/h**3/coef).tolist(),compensated_norm_max=float(co.max()),compensated_material_RMS_max=max(r['compensated_material_RMS'] for r in trials),compensated_log_slope='NOT_IDENTIFIABLE_FLOATING_FLOOR; do not fit order to constant numerical noise',operators=comp['operators'],interpretation='same-material gauge; cubic source-coordinate coefficient annihilated by two property directions; no new material branch')
 write('GAUGE_SLOPES.json',gauge);csvwrite('COMPENSATION_SOURCE.csv',[dict(**r,norm_over_h3=r['uncompensated_white_norm']/r['h']**3,coefficient_ratio=r['uncompensated_white_norm']/r['h']**3/coef) for r in trials])
 # Local prior-fibre mean shift relative to a likelihood-only Gaussian mode.
 localrows=[];fibre_ablation=[];fibre_arrays={}
 from numpy.polynomial.hermite import hermgauss
 for case,ref in refs.items():
  z=np.load(HERE/f'raw/{case}_material_GH21.npz');prob=z['weights']/z['weights'].sum();branch=prob@z['source_probability']
  v=z['V'];chol=np.linalg.cholesky(z['C']);zz=np.linalg.solve(chol,(v-z['mode']).T).T;gh,gw=hermgauss(21);ghweights=(gw[:,None]*gw[None,:]).ravel()/np.pi
  assert len(v)==441 and np.max(abs(zz-np.sqrt(2)*np.array([(a,b) for a in gh for b in gh])))<1e-9
  logq=-np.log(2*np.pi)-.5*np.linalg.slogdet(z['C'])[1]-.5*np.sum(zz**2,axis=1);logflat=np.log(ghweights)-z['phi']-logq;pf=np.exp(logflat-logflat.max());pf/=pf.sum();flatmean=pf@v;priormean=prob@v;sd=ref['material_ROI_SD']
  fibre_ablation.append(dict(observation=case,flat_material_density_mean_real=float(flatmean[0]),flat_material_density_mean_imag=float(flatmean[1]),mode_to_flat_likelihood_gap_SD=float(np.linalg.norm(z['mode']-flatmean)/sd),flat_likelihood_to_prior_fibre_gap_SD=float(np.linalg.norm(flatmean-priormean)/sd),mode_to_prior_fibre_gap_SD=ref['Gaussian_mode_mean_gap_SD'],scope='same saved GH21 nodes/responses, algebraic prior-density removal; local conditional Q1/Q2 only'))
  fibre_arrays[case+'_V']=v;fibre_arrays[case+'_flat_density_weights']=pf;fibre_arrays[case+'_prior_fibre_weights']=prob;fibre_arrays[case+'_flat_density_mean']=flatmean
  localrows.append(dict(observation=case,mean_real=ref['mean_increment'][0],mean_imag=ref['mean_increment'][1],mode_real=ref['material_Gaussian_mode'][0],mode_imag=ref['material_Gaussian_mode'][1],mode_mean_gap_abs=abs(cmean(ref)-complex(*ref['material_Gaussian_mode'])),mode_mean_gap_SD=ref['Gaussian_mode_mean_gap_SD'],material_SD=ref['material_ROI_SD'],host_mass=float(branch[0]),transition_mass=float(branch[1]),pure_mass=float(branch[2]),branch_label_disclaimer='quadrature probabilities for overlapping source strata, not distinct physical materials'))
 csvwrite('GN_FIBRE_REFERENCE.csv',localrows);csvwrite('PRIOR_FIBRE_REWEIGHTING.csv',fibre_ablation);np.savez_compressed(OUT/'PRIOR_FIBRE_REWEIGHTING.npz',**fibre_arrays)
 # Exact covariance and finite differences within the same saved finite node law.
 covrows=[];fdrows=[];savecov={}
 for method in ['S1_scaled_32_24','S3_source_16_9','S3_quotient_20_8']:
  z=np.load(HERE/f'raw/{method}.npz');V=z['V'];m=np.column_stack((V.real,V.imag));F=pack(z['prediction'])*np.sqrt(2);pw=z['prior_weights'];layer=z['layer'];pred=z['prediction']
  for oname,(y,sigma) in obs.items():
   phi=np.sum(abs(pred-y)**2,axis=(1,2))/sigma**2;w=pw*np.exp(-(phi-phi.min()));p=w/w.sum();mu=p@m;fm=p@F;C=(m-mu).T@(p[:,None]*(F-fm));between=np.zeros_like(C)
   for k in range(3):
    mask=layer==k;mass=p[mask].sum()
    if mass>1e-290:
     km=p[mask]@m[mask]/mass;kf=p[mask]@F[mask]/mass;between+=mass*np.outer(km-mu,kf-fm)
   within=C-between;D=C/sigma**2;covrows.append(dict(method=method,observation=oname,total_covariance_norm=float(np.linalg.norm(C)),within_covariance_norm=float(np.linalg.norm(within)),between_covariance_norm=float(np.linalg.norm(between)),total_derivative_norm=float(np.linalg.norm(D)),within_derivative_norm=float(np.linalg.norm(within)/sigma**2),between_derivative_norm=float(np.linalg.norm(between)/sigma**2),decomposition_residual=float(np.linalg.norm(C-within-between)),between_over_total_norm=float(np.linalg.norm(between)/max(np.linalg.norm(C),1e-300)),scope='finite saved-node conditional low-dimensional model; grid may be an unqualified posterior reference'))
   savecov[method+'_'+oname+'_total']=C;savecov[method+'_'+oname+'_between']=between
   if oname.startswith('HOST'):
    direction=np.linalg.svd(C,full_matrices=False)[2][0];score=F@direction;score-=score.mean();analytic=D@direction
    for relative in [1e-2,1e-3,1e-4]:
     delta=sigma*relative
     def tilted(sign):
      a=sign*delta*score/sigma**2;pp=p*np.exp(a-a.max());pp/=pp.sum();return pp@m
     fd=(tilted(1)-tilted(-1))/(2*delta);fdrows.append(dict(method=method,observation=oname,step_over_sigma=relative,analytic_directional_derivative_norm=float(np.linalg.norm(analytic)),finite_difference_norm=float(np.linalg.norm(fd)),relative_error=float(np.linalg.norm(fd-analytic)/max(np.linalg.norm(analytic),1e-30)),physics_RHS=0,qualification='same fixed response nodes reweighted, validates algebra only'))
 np.savez_compressed(OUT/'COVARIANCE_ARRAYS.npz',**savecov);csvwrite('WITHIN_BETWEEN_DERIVATIVES.csv',covrows);csvwrite('FIXED_NODE_FINITE_DIFFERENCES.csv',fdrows)
 # Receipt-grounded costs and qualification: no per-task recharging of shared nodes.
 globalreceipt=read(ROOT/'raw/jobs/stratified_uniform_submodels/CLOSED.json');localreceipt=read(ROOT/'raw/jobs/local_material_GN_reference/CLOSED.json');stages=[]
 for method in dict.fromkeys(r['method'] for r in rows):
  r=next(r for r in rows if r['method']==method);stages.append(dict(method=method,nodes=r['node_count'],incremental_physics_RHS=r['additional_shared_RHS'],counterfactual_uncached_forward_RHS=8*r['node_count'],first_summary_wall_seconds=r['method_wall_seconds'],scope='one response-node computation reused across all5observations; actual cache sharing includes host/material overlap'))
 csvwrite('GLOBAL_METHOD_COSTS.csv',stages)
 actual_receipt_RHS=sum(r['actual_counters']['actual_forward_RHS']+r['actual_counters']['actual_derivative_RHS'] for r in [globalreceipt,localreceipt]);assert actual_receipt_RHS==157656
 costs=dict(actual_A_total_RHS=actual_receipt_RHS,actual_receipt_forward_RHS=157304,actual_receipt_derivative_RHS=352,cap_RHS=160000,unused_RHS=2344,global_job_RHS=results['new_RHS'],local_job_RHS=157656-results['new_RHS'],global_job_wall_seconds=globalreceipt['wall_seconds'],local_job_wall_seconds=localreceipt['wall_seconds'],global_job_CPU_seconds=globalreceipt['CPU_seconds'],local_job_CPU_seconds=localreceipt['CPU_seconds'],local_GH15_node_RHS_each=225*8,local_GH21_node_RHS_each=441*8,local_GH15plus21_node_RHS_each=666*8,local_mode_J_line_search_RHS_total=(157656-results['new_RHS'])-2*666*8,task_cost_separation='NOT_ESTABLISHED',reason='host coarsest examined rules already pass both mean and mass refinement gates; regular global grids fail mean and overlapbranch gates; local reference has extra mode/J overhead and no global tail certificate',global_method_shared_incremental_RHS=sum(r['incremental_physics_RHS'] for r in stages),global_observation_setup_plus_compensation_RHS=results['new_RHS']-sum(r['incremental_physics_RHS'] for r in stages),new_analysis_physics_RHS=0,new_training_jobs=0)
 write('COST_AND_SCOPE.json',costs)
 # Figure1: host scaling and branch mixture.
 fig,axs=plt.subplots(1,3,figsize=(12,3.8))
 for method in ['S1_scaled_32_24','S3_source_16_9','S3_quotient_20_8']:
  rr=[r for r in host if r['method']==method];sx=[r['sigma_ratio'] for r in rr];label=method.replace('_',' ');axs[0].loglog(sx,[r['mean_magnitude'] for r in rr],'o-',label=label);axs[1].loglog(sx,[r['nonhost_mass'] for r in rr],'o-');axs[2].loglog(sx,[r['SD'] for r in rr],'o-')
 for ax,ylabel in zip(axs,['Conditional |mean V|','Nonhost source mass','Conditional material SD']):ax.set_xlabel('Noise / base noise');ax.set_ylabel(ylabel);ax.set_xticks([.25,.5,1.],['.25','.5','1']);ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(alpha=.2)
 axs[0].legend(fontsize=7);fig.suptitle('Host observations: same innovation at three noise scales; conditional submodels');fig.tight_layout();fig.savefig(OUT/'FIG01_HOST_SCALING.png',dpi=180);plt.close(fig)
 # Figure2: quadrature mean and label errors against independent local reference.
 fig,axs=plt.subplots(1,2,figsize=(11,4));methods=['S3_source_10_8','S3_source_16_9','S3_quotient_12_6','S3_quotient_20_8']
 for case,mark in [('REGULAR','o'),('PURE_OVERLAP','s')]:
  rr=[next(q for q in quals if q['comparison']==m+'->GN_GH21' and q['observation']==case) for m in methods];xx=np.arange(4);axs[0].semilogy(xx,[r['mean_error_SD'] for r in rr],mark+'-',label=case);axs[1].plot(xx,[r['branch_mass_max_abs_error'] for r in rr],mark+'-',label=case)
 axs[0].axhline(.01,color='black',ls=':',label='Frozen 1% SD gate');axs[1].axhline(.002,color='black',ls=':',label='Frozen .002 mass gate')
 for ax,ylabel in zip(axs,['Mean error / local reference SD','Max source-branch mass error']):ax.set_xticks(range(4),['Source10/8','Source16/9','Quotient12/6','Quotient20/8'],rotation=15,fontsize=8);ax.set_ylabel(ylabel);ax.legend(fontsize=8,loc='lower left');ax.grid(alpha=.2)
 fig.suptitle('Accurate prior normalization does not resolve narrow regular likelihoods');fig.tight_layout();fig.savefig(OUT/'FIG02_GLOBAL_GRID_MISS.png',dpi=180);plt.close(fig)
 # Figure3: gauge and exact cubic coefficient versus compensated floor.
 fig,axs=plt.subplots(1,2,figsize=(10,4));axs[0].loglog(h,un,'o-',label='Uncompensated response');axs[0].loglog(h,coef*h**3,'--',label='||c|| h^3');axs[0].loglog(h,co,'s-',label='Two-property compensation');axs[0].set_xticks([.005,.01,.02,.04,.08],['.005','.01','.02','.04','.08']);axs[0].xaxis.set_minor_formatter(NullFormatter());axs[0].set_xlabel('Distance h from upper seam');axs[0].set_ylabel('Whitened response change');axs[0].legend(fontsize=8);axs[0].grid(alpha=.2);axs[1].bar(['Real property only','Both properties'],[r['normal_ratio'] for r in comp['operators']]);axs[1].set_yscale('log');axs[1].set_ylabel('Projected cubic coefficient / ||c||');axs[1].set_title('Same-material gauge removes coefficient');fig.tight_layout();fig.savefig(OUT/'FIG03_GAUGE_COMPENSATION.png',dpi=180);plt.close(fig)
 # Figure4: resolved local conditional material density and source labels.
 fig,axs=plt.subplots(1,2,figsize=(11,4))
 for ax,case in zip(axs,['REGULAR','PURE_OVERLAP']):
  z=np.load(HERE/f'raw/{case}_material_GH21.npz');p=z['weights']/z['weights'].sum();v=z['V'];ref=refs[case];ax.scatter(v[:,0],v[:,1],c=np.log10(np.maximum(p,1e-30)),s=10+600*p/p.max(),cmap='viridis',alpha=.65);ax.plot(*ref['material_Gaussian_mode'],'rx',ms=10,label='Likelihood Gaussian mode');ax.plot(*fibre_arrays[case+'_flat_density_mean'],'bD',ms=5,label='Flat material-density mean');ax.plot(*ref['mean_increment'],'k+',ms=10,label='Prior-fibre corrected mean');ax.set_xlabel('Re V');ax.set_ylabel('Im V');ax.set_title(f"{case}: gap={ref['Gaussian_mode_mean_gap_SD']:.4g} SD");ax.legend(fontsize=7)
 fig.suptitle('Local GN-whitened GH21 quadrature; no global posterior/tail certificate');fig.tight_layout();fig.savefig(OUT/'FIG04_LOCAL_FIBRE_MEAN.png',dpi=180);plt.close(fig)
 # Save every figure's underlying sources explicitly in an archive manifest.
 write('FIGURE_SOURCES.json',{'FIG01_HOST_SCALING.png':['HOST_TRAJECTORIES.csv','HOST_SLOPES.json'],'FIG02_GLOBAL_GRID_MISS.png':['TASK_QUALIFICATION.csv','GLOBAL_GRID_TASKS.csv','GN_FIBRE_REFERENCE.csv'],'FIG03_GAUGE_COMPENSATION.png':['COMPENSATION_SOURCE.csv','GAUGE_SLOPES.json'],'FIG04_LOCAL_FIBRE_MEAN.png':['../raw/REGULAR_material_GH21.npz','../raw/PURE_OVERLAP_material_GH21.npz','GN_FIBRE_REFERENCE.csv','PRIOR_FIBRE_REWEIGHTING.csv','PRIOR_FIBRE_REWEIGHTING.npz']})
 summary=dict(status='OFFLINE_ANALYSIS_COMPLETE_Q1_Q2_ONLY',mean_and_branch_cost_separation='NOT_ESTABLISHED',global_tail_certificate=False,local_reference_refinement=[q for q in quals if q['method']=='GH15'],host_slopes=hostslopes,gauge=gauge,local_reference=localrows,prior_fibre_reweighting=fibre_ablation,costs=costs,max_fixed_node_FD_relative_error_at_step1e_4=max(r['relative_error'] for r in fdrows if r['step_over_sigma']==1e-4),physics_RHS=0,new_training_jobs=0,CPU_threads=1,analysis_wall_seconds=time.perf_counter()-tick,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 write('ANALYSIS_SUMMARY.json',summary)
 write('PLAIN_CLAIMS.json',[
  dict(claim='Host conditional mean decreases faster than nonhost source probability across the three fixed noise observations',qualification='Q1/Q2 conditional fixed2D/3D submodels, sameinnovation; no population rate proof',evidence=['HOST_TRAJECTORIES.csv','HOST_SLOPES.json']),
  dict(claim='Prior normalization near1 does not prevent global grids missing narrow regular likelihoods',qualification='allfour S3global rules fail frozen1%SD mean gate against localGN-fibre reference',evidence=['TASK_QUALIFICATION.csv','FIG02_GLOBAL_GRID_MISS.png']),
  dict(claim='Source transition/pure labels can overlap at the same material',qualification='local overlap GH21 masses .5688963452/.4311036548, no distinctmaterial-mode claim',evidence=['GN_FIBRE_REFERENCE.csv']),
  dict(claim='Prior fibre weighting shifts the regular local quadrature mean by .111978SD relative to flatmaterialdensity control',qualification='same savedGH21 finite nodes; not a full66 or globaltail result',evidence=['PRIOR_FIBRE_REWEIGHTING.csv','PRIOR_FIBRE_REWEIGHTING.npz']),
  dict(claim='Two adjustable properties annihilate the uniform material cubic source coefficient',qualification='exact same-material gauge and floating response floor; no new nontrivial material fibre',evidence=['GAUGE_SLOPES.json','COMPENSATION_SOURCE.csv']),
  dict(claim='Mean versus branch probability cost separation NOT_ESTABLISHED',qualification='coarsehost rules passboth; regularglobal meanfails; localrefmode/Joverheadcharged; no globaltailcert',evidence=['COST_AND_SCOPE.json','TASK_QUALIFICATION.csv']),
  dict(claim='Within/between covariance identity matches saved-node directional finite differences',qualification='algebra check within samefinite quadrature rule only, including potentially unqualified globalposteriorrules',evidence=['WITHIN_BETWEEN_DERIVATIVES.csv','FIXED_NODE_FINITE_DIFFERENCES.csv'])])
 # Narrative scientific claims remain parent-reviewed.
 text=['# Track A：条件子模型的材料商、任务积分与实际成本','', '本分析只读取已经付费并保存的response/积分节点和receipt，新增物理RHS=0，新增训练=0。实际A共157656/160000 RHS；全局网格146552、局部GN白化参考11104。研究对象是固定2D/3D、active Gaussian截断[-4,4]的条件子模型，**不是Full66边缘后验**。','', '## 主要可复核结果','', 'host三个固定噪声尺度下，材料均值随噪声缩小，host源质量增加。S3 source finest的σ比例1/.5/.25对应的均值实部为8.06050e-6、3.18995e-6、1.26501e-6，nonhost质量为.0590994、.0471467、.0375571。三点受控斜率只说明这组三个条件观测，不是总体Bayes风险率。','', '全局网格的prior normalization准确，不保证窄likelihood被解析。四个S3 global规则对REGULAR/PURE_OVERLAP均未通过固定1%SD材料均值门。PURE_OVERLAP时它们几乎全分配给transition源标签，但独立局部GN+prior-fibre GH21积分为transition=.5688963452、pure=.4311036548；两个源标签可对应同一物理材料，不能称为两个已区分的材料世界。','', '局部GN白化GH15→21的REGULAR均值差3.93e-10、overlap约1e-15；材料均值、branch、latent/material事件均通过已冻结的1e-5绝对/1%SD/.002 refinement门。它们仍只有Q1/Q2数值收敛证据：高斯重要性网格的尾部和外部其他mode没有全局上界。没有定义或声称离散节点测度对连续密度的TV精度。','', 'Likelihood-only Gaussian mode与保留prior fibre volume的积分均值之间，REGULAR相差.10936个材料SD，overlap相差.002647个SD。mode与积分均值的比较本身混合了似然非线性和prior fibre权重，不能将likelihoodmode等同于完整Laplace posterior均值。另对同一批已存GH21节点移除prior材料密度（没有新查询），REGULAR mode→flat-likelihood mean仅.002626SD，而flat-likelihood→prior-fibre mean为.111978SD；overlap对应.002565SD与.000617SD。该固定节点对照更直接显示REGULAR中的prior fibre volume作用，仍仅是局部条件模型结果。','', '上接缝未补偿response的log斜率接近3；两项属性补偿后的材料变化停留于1.48e-16 RMS、白化response约4.15e-13浮点底噪。对底噪拟合高阶斜率没有意义。只开放实属性时，投影后仍保留8.0983%的立方系数；实虚属性均开放时降至1.56e-15相对比例。这是同材料gauge消去源坐标系数，没有建立新的非平凡材料纤维。','', '## 冻结任务门与成本解释','', '逐任务门保留原值：材料均值绝对差≤1e-5且差/参考材料SD≤.01；source branch、latent phase、material事件概率差各≤.002。SD自身与logevidence差另报，不把它们换成更容易的成功标准。所有状态、失败和确切数值见TASK_QUALIFICATION.csv。','', 'host粗网格在当前refinement比较中已同时满足均值与mass门，regular全局规则没有达到材料均值门；局部参考又增加mode/J/line-search开销。因此**材料均值比branch概率更便宜的成本分离没有建立**，不能把共享一次节点的多任务结果重复收费或从单一host示例推广。全局每规则actual增量RHS与uncached反事实分别保存；本轮没有扩展到4D以救结论。','', 'within/between task-response covariance及observation导数使用现有同一节点规则重新加权计算。导数坐标为sqrt(2)×source-major pack(y)，所以导数矩阵为Cov(V,sqrt(2)pack(F))/sigma²。host层内材料固定仍可因host/transition权重迁移产生between项。中央差分只重新评价存储节点的likelihood，新增RHS为0；它验证有限规则的代数导数，不补足该规则的posterior reference资格。','', '## 文件与图','', 'ANALYSIS_SUMMARY.json、TASK_QUALIFICATION.csv、GLOBAL_GRID_TASKS.csv、HOST_TRAJECTORIES.csv、WITHIN_BETWEEN_DERIVATIVES.csv、FIXED_NODE_FINITE_DIFFERENCES.csv、COVARIANCE_ARRAYS.npz、COST_AND_SCOPE.json及FIGURE_SOURCES.json给出所有数据和口径。','', '![host条件尺度](FIG01_HOST_SCALING.png)','','![全局网格遗漏窄likelihood](FIG02_GLOBAL_GRID_MISS.png)','','![gauge补偿](FIG03_GAUGE_COMPENSATION.png)','','![局部prior-fibre均值](FIG04_LOCAL_FIBRE_MEAN.png)']
 (OUT/'REPORT_ZH.md').write_text('\n'.join(text)+'\n')
 print(json.dumps(dict(status=summary['status'],analysis_seconds=summary['analysis_wall_seconds'],host_slopes=hostslopes,gauge_slope=gauge['uncompensated_last3_log_slope'],local_reference=localrows,FD_max=summary['max_fixed_node_FD_relative_error_at_step1e_4']),indent=2))
if __name__=='__main__':main()
