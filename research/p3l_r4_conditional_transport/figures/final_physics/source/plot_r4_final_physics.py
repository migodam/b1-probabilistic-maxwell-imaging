#!/usr/bin/env python3
"""Closed-array physical evidence panels. No Target/backend/GT/model creation."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write
from p3r4.evaluation import checkpoint_identity
DEST=ROOT/'figures/final_physics'
def main():
 if DEST.exists():raise FileExistsError('preserve first figure package')
 b=Book('IV','final_physics_figures_v1',max_wall=900,analysis_only=True);status='COMPLETE';error=None
 try:
  import numpy as np
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  from matplotlib.ticker import FuncFormatter
  plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','axes.labelcolor':'#26374a','text.color':'#26374a','axes.edgecolor':'#b1bac5','pdf.fonttype':42,'savefig.facecolor':'white'})
  DEST.mkdir(parents=True);sources={}
  def load(p):
   p=Path(p);sources[str(p)]=checkpoint_identity(p);return json.loads(p.read_text())
  def closed(p):
   d=load(p)
   if d['status']!='COMPLETE':raise ValueError('closed only:'+str(p))
   return d
  def table(name,rows):
   with (DEST/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  def save(fig,name):
   fig.savefig(DEST/f'{name}.png',dpi=220);fig.savefig(DEST/f'{name}.pdf');plt.close(fig)
  C=['#526579','#2b867e','#538fb4','#b68b42'];arms=('full','richardson','global_gmres','global_gmres_accuracy');labels=('Full','Richardson','Global GMRES','GMRES fixed accuracy')
  base=OUT/'results/IV/interval_chain_locked_v1';comparison=load(base/'COMPARISON.json');closed(OUT/'results/IV/interval_chain_locked_pair/COMPLETE.json')
  if comparison['total_misses']!=0:raise ValueError('source comparison inconsistent')
  chain=[];chainparents=[]
  for arm,label in zip(arms,labels):
   receipt=closed(OUT/f'IV/jobs/interval_chain_locked_v1_{arm}/CLOSED.json');ctr=receipt['counters'];setup=0;accept=fall=0;setupwall=0
   for i in range(8):
    s=load(base/arm/f'R4L1_{i:05d}/SUMMARY.json');accept+=s['accepts'];fall+=s['fallbacks'];ss=s.get('setup',{});setup+=ss.get('RHS',0);setupwall+=ss.get('wall_including_R0_boundcopy_and_persistence',0)
    chainparents.append(dict(parent=s['parent'],arm=arm,steps=s['steps'],accepts=s['accepts'],fallbacks=s['fallbacks'],RHS=s['costs']['RHS'],setup_RHS=ss.get('RHS',0),wholeparentwall=s['wholeparentwall']))
   cost=comparison['costs'][arm];small=ctr.get('reduced_LS_RHS',0);forward=ctr.get('full_forward_RHS',0)
   assert cost['actual_RHS']==setup+small+forward
   chain.append(dict(arm=arm,label=label,all_complex_RHS=cost['actual_RHS'],full_forward_RHS=forward,R0_setup_RHS=setup,reduced_LS_RHS=small,whole_process_seconds=cost['whole_process_seconds'],controlled_seconds=receipt['controlled_wall_seconds'],CPU_seconds=receipt['profile']['CPU_seconds'],R0_setup_wall_seconds=setupwall,steps=8192,accepts=accept,fallbacks=fall,new_states=ctr.get('full_forward_calls',0),new_LU=ctr.get('LU_factorizations',0)))
  table('chain_costs.csv',chain);table('chain_per_parent.csv',chainparents)
  fig,ax=plt.subplots(1,2,figsize=(11.7,4.9));fig.subplots_adjust(left=.08,right=.98,bottom=.24,top=.82,wspace=.32);fig.suptitle('Same-cache fixed chain: measured cost on eight new locked parents',y=.96,fontsize=15,fontweight='bold')
  yy=np.arange(4);vals=[r['whole_process_seconds'] for r in chain]
  ax[0].barh(yy,vals,color=C,height=.58);ax[0].scatter([r['controlled_seconds'] for r in chain],yy,marker='|',s=160,color='white',linewidth=2,label='Book controlled wall')
  for i,v in enumerate(vals):ax[0].text(v+1,i,f'{v:.2f} s',va='center',fontsize=10)
  ax[0].set(yticks=yy,yticklabels=labels,xlabel='Whole arm process wall (seconds)',xlim=(0,85));ax[0].invert_yaxis();ax[0].set_title('A  Includes I/O, setup and arm postprocessing',loc='left',fontsize=11);ax[0].legend(loc='lower right',fontsize=8,frameon=False)
  left=np.zeros(4)
  for key,label,col in [('full_forward_RHS','Full forward','#526579'),('R0_setup_RHS','R0 setup','#2b867e'),('reduced_LS_RHS','Reduced complex least squares','#b68b42')]:
   val=np.array([r[key] for r in chain]);ax[1].barh(yy,val,left=left,height=.58,label=label,color=col);left+=val
  for i,v in enumerate(left):ax[1].text(v+900,i,f'{int(v):,}',va='center')
  ax[1].set(yticks=yy,yticklabels=[],xlabel='All charged complex right-hand sides',xlim=(0,80000));ax[1].invert_yaxis();ax[1].xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}k'));ax[1].set_title('B  Reduced-system RHS are included',loc='left',fontsize=11);ax[1].legend(loc='lower right',fontsize=8,frameon=False)
  fig.text(.08,.105,'8 parents × 1,024 pCN transitions; beta = .002. Same proposals/uniforms: zero saved path, decision or proof misses.',fontsize=9)
  fig.text(.08,.055,'Cached anchor/receiver rights shared by all arms; R0 setup is paid in each approximate arm. One quiet serial run; no ESS or posterior claim.',fontsize=8.5);save(fig,'fig01_locked_chain_cost')
  # Geometry: all planned attempts retained, including incomplete paths.
  geo=load(OUT/'results/offline_only/geometry_matched_gate_audit_v1/RESULTS.json');proto=load(OUT/'results/offline_only/geometry_matched_gate_audit_v1/FROZEN_ANALYSIS_PROTOCOL.json')
  for job in ('dynamic_geometry_first','dynamic_geometry_heun_r1'):closed(OUT/f'III/jobs/{job}/CLOSED.json')
  closed(OUT/'III/jobs/geometry_matched_gate_audit_v1/CLOSED.json')
  order=[('dynamic_geometry_first','direct'),('dynamic_geometry_first','frozen'),('dynamic_geometry_first','dynamic'),('dynamic_geometry_first','Euler'),('dynamic_geometry_heun_r1','direct'),('dynamic_geometry_heun_r1','Heun')]
  gt={ (r['campaign'],r['family']):r for r in geo['totals']};grows=[gt[k] for k in order];gnames=['Direct bracket (first)','Frozen weak (first)','Dynamic weak (first)','Four Euler segments','Direct bracket (revision)','Four Heun segments']
  table('geometry_total_gates.csv',[dict(label=n,**r) for n,r in zip(gnames,grows)]);table('geometry_rank_radius_dimension.csv',geo['groups'])
  fig=plt.figure(figsize=(12.4,8.6));gs=fig.add_gridspec(2,1,height_ratios=[1.1,1.1],left=.22,right=.95,bottom=.14,top=.86,hspace=.66);a=fig.add_subplot(gs[0]);hax=fig.add_subplot(gs[1]);fig.suptitle('Finite material-distance gates: keep every planned attempt',y=.96,fontsize=15,fontweight='bold')
  v=np.array([[r['all_support_and_matched'],r['matched_2pct']-r['all_support_and_matched'],r['completed']-r['matched_2pct'],r['planned']-r['completed']] for r in grows]);den=np.array([r['planned'] for r in grows]);left=np.zeros(6)
  for j,(lab,col) in enumerate(zip(['Matched + all support','Matched, support incomplete','Completed but unmatched','Incomplete'],['#2b867e','#adcac6','#d5dce4','#b68b42'])):
   x=v[:,j]/den;a.barh(range(6),x,left=left,color=col,height=.66,label=lab);left+=x
  a.set(yticks=range(6),yticklabels=gnames,xlim=(0,1.30),xticks=[0,.25,.5,.75,1],xlabel='Fraction of all planned attempts');a.invert_yaxis();a.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:.0%}'))
  for i,r in enumerate(grows):a.text(1.025,i,f"{r['matched_2pct']}/{r['planned']} matched\n{r['all_support_and_matched']}/{r['planned']} all support",va='center',fontsize=8.5)
  a.set_title('A  Exact counts use the planned denominator',loc='left',fontsize=11);a.legend(loc='upper center',bbox_to_anchor=(.42,-.28),ncol=2,fontsize=8,frameon=False)
  pairs=[(rank,radius) for rank in (2,4,8) for radius in (.02,.04,.08)];grid=np.zeros((6,9));texts=[]
  for i,k in enumerate(order):
   tex=[]
   for j,(rank,radius) in enumerate(pairs):
    sel=[r for r in geo['groups'] if (r['campaign'],r['family'])==k and r['rank']==rank and r['radius']==radius];n=sum(r['planned_attempts'] for r in sel);m=sum(r['all_support_and_match_2pct_count'] for r in sel);grid[i,j]=m/n;tex.append(f'{m}/{n}')
   texts.append(tex)
  im=hax.imshow(grid,cmap='YlGnBu',vmin=0,vmax=1,aspect='auto');hax.set(yticks=range(6),yticklabels=gnames,xticks=range(9),xticklabels=[f'k={k}\nr={r:.02f}' for k,r in pairs]);hax.tick_params(axis='both',length=0)
  for i in range(6):
   for j in range(9):hax.text(j,i,texts[i][j],ha='center',va='center',fontsize=8.5,color='white' if grid[i,j]>.55 else '#26374a')
  hax.set_title('B  Matched + all support by rank and target RMS (64/66 dimensions pooled)',loc='left',fontsize=11,pad=12)
  fig.text(.035,.065,'8 exposed DEV parents; repeated settings and orientations are attempts, not independent worlds. Distance gate: |actual RMS / target − 1| ≤ .02.',fontsize=8.5)
  fig.text(.035,.035,'Support = original noise ball, ΔPhi ≤ 2, prior Δ ≤ 2, signed progress ≥ .5r, prior legality and residual check. Sampled paths give no continuous certificate.',fontsize=8.2);save(fig,'fig02_geometry_gates')
  # Risk/data-fit tradeoff: reuse already scored fixed64 loss arrays and paid Full Phi.
  scorebase=OUT/'results/offline_only/L1_optimization_addenda_GT_v1';score=load(scorebase/'RESULTS.json');closed(OUT/'IV/jobs/L1_optimization_addenda_isolated_GT_v1/CLOSED.json')
  scorefile=scorebase/'PER_PARENT_LOSSES.npz';sources[str(scorefile)]=checkpoint_identity(scorefile)
  with np.load(scorefile,allow_pickle=False) as z:actions=z['actions'].tolist();mse=z['full_mse'].copy();parents=z['parent'].tolist()
  raw=load(OUT/'results/fresh/locked/raw_forward/RESULTS.json');closed(OUT/'IV/jobs/fresh_raw_locked_v1/CLOSED.json');closed(OUT/'IV/jobs/fresh_locked64_representatives_v1/CLOSED.json')
  rd={(r['parent'],r['action']):r for r in raw['rows']};pp=[];means=[];keys=('spline_raw','spline_clipped','r2_raw','r2_clipped')
  for i,p in enumerate(parents):
   anchorphi=None
   for key in keys:
    rep=load(OUT/f'results/fresh/locked/representatives/{p}/{key}/RESULTS.json');anchorphi=rep['anchor_Phi']
    for kind,act,phi in [('raw_risk','fixed_'+key,rd[(p,key)]['Phi']),('representative','original_rep_'+key,rep['Phi'])]:
     pp.append(dict(parent=p,key=key,kind=kind,Phi=phi,full_mse=mse[i,actions.index(act)],data_feasible=phi<=282.873792,source_role='rawdecoderrealizabilityNOT_PROVEN' if kind=='raw_risk' else 'finite_budget_candidate_not_posterior'))
   pp.append(dict(parent=p,key='anchor',kind='anchor',Phi=anchorphi,full_mse=mse[i,actions.index('original_anchor')],data_feasible=anchorphi<=282.873792,source_role='online_anchor'))
  for key in (*keys,'anchor'):
   for kind in (('anchor',) if key=='anchor' else ('raw_risk','representative')):
    ss=[r for r in pp if r['key']==key and r['kind']==kind];means.append(dict(key=key,kind=kind,n_parent=len(ss),mean_Phi=float(np.mean([r['Phi'] for r in ss])),mean_full_mse=float(np.mean([r['full_mse'] for r in ss])),noiseball_count=sum(r['data_feasible'] for r in ss)))
  table('risk_fit_per_parent.csv',pp);table('risk_fit_means.csv',means)
  fig,ax=plt.subplots(1,2,figsize=(11.6,5.2));fig.subplots_adjust(left=.085,right=.97,bottom=.23,top=.82,wspace=.28);fig.suptitle('Material risk and nonlinear data fit are different quantities',y=.96,fontsize=15,fontweight='bold')
  anchor=means[-1]
  for j,clip in enumerate(('raw','clipped')):
   a=ax[j];a.axvline(282.873792,color='#a97842',linestyle='--',linewidth=1,label='Original noise-ball threshold')
   for fam,col,dy in [('spline','#2b867e',.0005),('r2','#538fb4',-.00045)]:
    key=fam+'_'+clip;rr=next(r for r in means if r['key']==key and r['kind']=='raw_risk');re=next(r for r in means if r['key']==key and r['kind']=='representative')
    a.annotate('',xy=(re['mean_Phi'],re['mean_full_mse']),xytext=(rr['mean_Phi'],rr['mean_full_mse']),arrowprops=dict(arrowstyle='->',color=col,lw=1.6))
    a.scatter(rr['mean_Phi'],rr['mean_full_mse'],s=76,marker='s',color=col,zorder=3);a.scatter(re['mean_Phi'],re['mean_full_mse'],s=76,marker='o',facecolor='white',edgecolor=col,lw=2,zorder=3)
    a.annotate(f"{fam.title()} risk ({rr['noiseball_count']}/64 in ball)",(rr['mean_Phi'],rr['mean_full_mse']),xytext=(3,8 if fam=='spline' else -17),textcoords='offset points',fontsize=8,color=col)
    a.annotate(f"{fam.title()} representative",(re['mean_Phi'],re['mean_full_mse']),xytext=(7,9 if fam=='spline' else -16),textcoords='offset points',fontsize=8,color=col)
   a.scatter(anchor['mean_Phi'],anchor['mean_full_mse'],s=70,color='#526579',marker='D',zorder=4);a.annotate('Original anchor',(anchor['mean_Phi'],anchor['mean_full_mse']),xytext=(7,4),textcoords='offset points',fontsize=8)
   a.set_title(('A  Raw estimators','B  Component-clipped estimators')[j],loc='left',fontsize=11);a.set_xlabel('Mean nonlinear Full Phi');a.set_ylabel('Mean full 576-cell material MSE');a.grid(alpha=.13);a.set_ylim(.0056,.0114);a.margins(x=.30)
  ax[0].legend(loc='lower left',fontsize=8,frameon=False)
  fig.text(.085,.115,'Identical first 64 L1 parents; fixed Spline seed 2026102501 and frozen R2. Squares: material risk actions; circles: original bounded representatives.',fontsize=8.2)
  fig.text(.085,.075,'Arrows link fixed outputs, not an optimization trajectory. Raw/clipped Full evaluation does not prove decoder realizability; representatives are not posterior means.',fontsize=8.2)
  fig.text(.085,.035,'All 320 original SLSQP cases stopped at budget. Callback-only selection limitation is preserved; independent Full-J controls are shown separately.',fontsize=8.2);save(fig,'fig03_risk_fit_tradeoff')
  # Exact Full-J diagnostics from already-closed TR controls, no derivative recomputation.
  trbase=OUT/'results/fresh/locked/ordinary_fullJ_TR_L1_64_addendum_v1';tr=load(trbase/'RESULTS.json');closed(OUT/'IV/jobs/ordinary_fullJ_TR_L1_64_addendum_v1/CLOSED.json');diag=[]
  for case in tr['cases']:
   diag.append(dict(parent=case['parent'],lambda_=case['lambda_'],initial_GN=case['initial_kkt']['GN_Newton_decrement'],final_GN=case['kkt']['GN_Newton_decrement'],initial_gradient_inf=case['initial_kkt']['gradient_total_inf'],final_gradient_inf=case['kkt']['gradient_total_inf'],delta_Phi=case['Phi']-case['initial_kkt']['Phi'],delta_objective=case['objective']-case['initial_kkt']['objective'],RHS=case['actual_cost']['RHS'],accepted_updates=case['accepted_updates'],rejected_trials=case['rejected_trials']))
  table('fullJ_TR_diagnostics.csv',diag)
  fig,ax=plt.subplots(1,2,figsize=(11.6,5.5));fig.subplots_adjust(left=.09,right=.96,bottom=.23,top=.82,wspace=.34);fig.suptitle('Exact Full-J controls: objective progress and GN stationarity differ',y=.96,fontsize=15,fontweight='bold')
  for lam,col,offset in [(0,'#538fb4',0),(1,'#2b867e',2.4)]:
   ss=[r for r in diag if r['lambda_']==lam];x=np.array([offset,offset+1])
   for r in ss:ax[0].plot(x,[r['initial_GN'],r['final_GN']],color=col,alpha=.16,lw=.7)
   med=[np.median([r[k] for r in ss]) for k in ('initial_GN','final_GN')];ax[0].plot(x,med,'o-',color=col,lw=2.4,label=f'lambda {lam} median')
   for xx,v in zip(x,med):ax[0].annotate(f'{v:.3g}',(xx,v),xytext=(0,9),textcoords='offset points',ha='center',fontsize=9,color=col)
   for metric,dx,m in [('delta_Phi',-.13,'o'),('delta_objective',.13,'x')]:
    yy=np.array([r[metric] for r in ss]);jit=np.linspace(-.055,.055,len(yy));ax[1].scatter(np.full(len(yy),lam)+dx+jit,yy,s=17,color=col,marker=m,alpha=.5,label=metric if lam==0 else None)
  ax[0].set(yscale='log',xticks=[0,1,2.4,3.4],xticklabels=['Initial\nλ=0','Final\nλ=0','Initial\nλ=1','Final\nλ=1'],ylabel='GN Newton decrement (log scale)',xlim=(-.4,3.8));ax[0].set_title('A  Every parent retained; lines are paired',loc='left',fontsize=11);ax[0].grid(axis='y',alpha=.15)
  ax[1].axhline(0,color='#929ca8',lw=1);ax[1].set(xticks=[0,1],xticklabels=['λ=0: data Phi','λ=1: Phi + prior'],ylabel='Final minus initial objective / Phi',yscale='symlog');ax[1].set_yscale('symlog',linthresh=1e-5);ax[1].legend(fontsize=8,frameon=False,loc='lower right');ax[1].set_title('B  Full nonlinear objective changes',loc='left',fontsize=11);ax[1].grid(axis='y',alpha=.15)
  fig.text(.085,.115,'64 paired L1 parents; post-method-freeze, pre-GT exploratory addendum. Exact Full-J gradient with GN curvature, not the full nonlinear Hessian.',fontsize=8.3)
  fig.text(.085,.075,'All 128 arms stopped at the same 16-state / 256-RHS local cap. Total new 30,752 RHS = Full 9,184 + accepted-update J 21,568; cached initial J disclosed.',fontsize=8.3)
  fig.text(.085,.035,'Small median decrement is not universal or global convergence: λ=1 final maximum = .542; exact gradient maximum = 2.548. No MMSE/posterior certificate.',fontsize=8.3);save(fig,'fig04_fullJ_TR_diagnostics')
  # Freeze source identities and exact scope, then wait for a visual QA marker.
  first8context=sum(load(OUT/f'results/contexts/locked/R4L1_{i:05d}.json')['new_RHS'] for i in range(8))
  facts=dict(status='FIGURES_WRITTEN_AWAIT_VISUAL_QA',new_RHS=0,new_Maxwell=0,new_LU=0,new_fit=0,new_GT_reads=0,source_identities=list(sources.values()),chain=chain,chain_common_first8_context_owner_RHS=first8context,geometry=grows,geometry_protocol=proto,risk_fit=means,scope='closedarraydescriptiveonly,noaddedscientificverdict',figure_files=sorted(p.name for p in DEST.glob('fig*')),old_source_bytes_preserved=True)
  write(DEST/'FACTS.json',facts)
  draft=f'''# 物理证据图说明草稿

这些图只读取封账数据与已封存评分；不新增Maxwell、模型、拟合或真值读取。英文图内文字方便论文排版，以下说明保留适用范围。

**图1：同缓存、同链的实际费用。** 新锁定前8对象，每对象1024步pCN、beta=.002；四臂使用相同随机数和缓存权限，源复核为零决策/路径/proof差异。左图使用完整进程墙钟，白线是Book受控墙钟；它们都不等于纯物理时间。右图逐项保留Full、R0建立和降阶复数最小二乘RHS。四臂总计105594 RHS。共同上下文前8对象原拥有{first8context} RHS，均另有owner，不按新图重复收费；每个近似臂另付R0的4608 RHS。此图是一轮安静串行固定链的数值和费用描述，不是ESS增益或后验参考结论。

**图2：有限材料距离与支持门槛。** 前8个已暴露DEV对象，计划尝试分母含所有方向、64/66维、秩2/4/8、半径.02/.04/.08，重复设置不等于独立世界。顶部明确分别显示距离合格及全部支持合格；下图每格保留成功数/计划数。2%距离门槛是封账后统一诊断，不改原算法/teacher 1%门槛；首轮套用修订版支持门槛仅描述，不追认原教师资格。直接法两campaign分开，不能当两个独立复现。四Euler自身场二阶缺陷不作为交换子理论反证；Heun失败也不作GT距离不匹配比较。noiseball用原几何协议，priorΔ≤2、ΔPhi≤2、signed progress≥.5r、合法性和残差检查均保留；采样路径不提供连续证书。

**图3：材料风险与数据拟合的不同角色。** 相同L1前64对象，固定Spline_s1及R2，直接读取已完成隔离评分的材料MSE与既有Full Phi；没有新GT访问。点是对象均值，不是根据GT选择的案例。raw/clipped风险动作虽经Full诊断，decoder可实现性仍NOT_PROVEN；圆点为原预算受限的代表元。箭头只连接两个固定输出，不能当优化轨迹。原320例均预算停止，callback返回选择的语义局限不隐藏。原查询池80640 RHS及L1 raw完整诊断24576 RHS归原receipt，画图不重计；此图的64子集不是完整512对象主评分。

**图4：Full-J强控制的诊断范围。** 仅展示GT前封存的λ0/λ1固定64对象附加控制。每个接受更新刷新J32，拒绝Full8也收费，新30752 RHS；初始J来自已付共同上下文。λ0的Phi下降和λ1的小GN decrement必须区分；GN矩阵不等于完整非线性Hessian，最大梯度/decrement及全部预算停止仍保留。此补充发生于方法冻结后、GT前，属于探索性正确度诊断，不改primary输出，也不证明全局MAP、MMSE或后验均值。

PNG/PDF、compact CSV、完整来源路径/bytes/mtime记录在本目录。所有原始源文件保持原字节；最终科学解释由主线程负责。
'''
  (DEST/'DRAFT_ZH.md').write_text(draft);print('FIGURES_READY '+str(DEST),flush=True)
  deadline=time.monotonic()+600
  while not (DEST/'VISUAL_QA.json').exists():
   if time.monotonic()>deadline:raise TimeoutError('visual QA not acknowledged')
   time.sleep(.5)
  qa=json.loads((DEST/'VISUAL_QA.json').read_text())
  if qa['status']!='PASS':raise ValueError('visual QA did not pass')
  for ident in sources.values():
   if checkpoint_identity(ident['path'])!=ident:raise ValueError('source changed:'+ident['path'])
  facts['status']='COMPLETE_VISUALLY_INSPECTED';facts['visual_QA']=qa;write(DEST/'FACTS.json',facts)
  print('COMPLETE zero_new_physics',flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
