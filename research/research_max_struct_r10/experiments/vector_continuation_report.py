"""Public-only incremental budget curves/report; never reads inverse GT."""
import time
T=time.perf_counter()
import json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'RESULTS/vector_continuation_v04';figs=ROOT/'FIGURES/vector';figs.mkdir(exist_ok=True)
Y=np.load(ROOT/'data_public/vector/v04_matched12/observations.npz',allow_pickle=False)['Y'];rows=[]
fig,ax=plt.subplots(figsize=(6.4,4))
for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
 p=out/arm;d=json.loads((p/'RECEIPT.json').read_text());c=d['physical_corrector'];old=ROOT/'RESULTS/vector'/('v04_'+arm);start=json.loads((old/'physical_accepted80.json').read_text());trace=[json.loads(x) for x in (p/'TRACE.jsonl').read_text().splitlines()];accepted=[x for x in trace if x.get('kind')=='physical_accepted']
 z=np.load(p/'prediction.npz',allow_pickle=False);global_res=float(np.linalg.norm(z['prediction']-Y)/np.linalg.norm(Y));z.close()
 ax.plot([0]+[x['accepted'] for x in accepted],[start['Phi']]+[x['Phi'] for x in accepted],label=arm)
 rows.append(dict(arm=arm,prefix_m80_Phi3=start['Phi'],continued_Phi3=d['physical_Phi'],source_relative_RMS_residual=float(np.sqrt(d['physical_Phi'])),global_Frobenius_relative_data_residual=global_res,PG_inf=d['PG_inf'],accepted=c['accepted'],requests=c['requests'],stop=c['stop'],nested_CPU_seconds=d['CPU_seconds'],nested_wall_seconds=d['wall_seconds'],physical_LU=d['counts']['physical_LU'],primal_RHS=d['counts']['primal_RHS'],adjoint_RHS=d['counts']['adjoint_RHS']))
ax.set(xlabel='Additional accepted physical steps (history reset)',ylabel='Public source-normalized Phi3',title='v04 posthoc budget diagnostic: same physical objective')
ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(alpha=.25);ax.legend();fig.tight_layout();figure=figs/'V04_CONTINUATION_PUBLIC_PHI.png';fig.savefig(figure,dpi=170);plt.close(fig)
with (out/'PUBLIC_ENDPOINTS.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
analysis=dict(platform='Mac',CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-T,scope='Whole independent public-only endpoint residual analysis and plot process; zero scientific RHS or GT reads',physical_RHS=0,GT_read=False)
(out/'REPORT_ANALYSIS_COST.json').write_text(json.dumps(analysis,indent=2)+'\n')
p=out/'ACCOUNTING_FINAL.json';ledger=json.loads(p.read_text());ledger['rows']=[x for x in ledger['rows'] if x['identity']!='public-report-analysis'];ledger['rows'].append(dict(identity='public-report-analysis',platform='Mac',CPU_seconds=analysis['CPU_seconds'],wall_seconds=analysis['wall_seconds'],source='RESULTS/vector_continuation_v04/REPORT_ANALYSIS_COST.json',scope=analysis['scope'],counts=dict(physical_RHS=0),additive=True,included_in=None));ledger['by_platform_CPU_lower_bounds']={host:sum(x['CPU_seconds'] for x in ledger['rows'] if x['platform']==host and x['CPU_seconds'] is not None) for host in ('Mac','XINAN')};p.write_text(json.dumps(ledger,indent=2)+'\n')
lines=['# v04 追加原物理预算诊断','', '分类：POSTHOC_BUDGET_DIAGNOSTIC，已曝光development，不能充作独立确认。三臂从各自原physical_accepted80公开材料出发，一致重置L-BFGS历史；原18预测和canonical不变，未重新造Y/噪声或重跑aux。工人未读GT。','', '目标保持Phi3=mean_s ||S J_s−Y_s||²/||Y_s||²，无1/2、R0、χbox[0,8]；同一G/E/S、60源。它不是全Y归一化，也不是iid同σ白化ML。','', '|arm|m80 Phi3|追加Phi3|PG∞|追加accepted/requests|global相对数据残差|','|---|---:|---:|---:|---:|---:|']
for d in rows:lines.append(f"|{d['arm']}|{d['prefix_m80_Phi3']:.9g}|{d['continued_Phi3']:.9g}|{d['PG_inf']:.3g}|{d['accepted']}/{d['requests']}|{d['global_Frobenius_relative_data_residual']:.6g}|")
lines+=['','三臂都以新增80步预算上限终止，PG仍约5.4–9.1e−6，不能称收敛/驻点或唯一性。公开Phi继续下降只说明原80步不是最终物理最小值，不能据此推断材料恢复；材料评分由root在三预测共同seal之后独立做。选择追加v04的依据是公开PG/fit，不是GT排名。','',f"新增已测CPU下界：Mac {ledger['by_platform_CPU_lower_bounds']['Mac']:.6f}s，XINAN {ledger['by_platform_CPU_lower_bounds']['XINAN']:.6f}s。一整GPU科学子进程含coldsetup/三臂试探/保存，计243LU+14580正向/14580伴随RHS；adjoint复用每个材料的原LU。query子CPU和transport另计，nested三臂/设置/stream-event不重加；原18prefix为非additive引用。手工实现与未单计progress/query等仍null。",'', '图：FIGURES/vector/V04_CONTINUATION_PUBLIC_PHI.png；费用：RESULTS/vector_continuation_v04/ACCOUNTING_FINAL.json；公开表：同目录PUBLIC_ENDPOINTS.csv。']
(ROOT/'delegated/vector/V04_CONTINUATION_REPORT_ZH.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status='COMPLETE_PUBLIC_ONLY',rows=rows,by_platform=ledger['by_platform_CPU_lower_bounds'])))
