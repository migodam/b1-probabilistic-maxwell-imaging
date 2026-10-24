#!/usr/bin/env python3
"""Read-only float32 boundary representation audit; no material alteration."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,save
from p1_material.prior import MaterialPrior
from p3r4.evaluation import checkpoint_identity

def main():
 dest=OUT/'AUDIT/raw_risk_actions_DEV_float32'
 if dest.exists():raise FileExistsError('preserve existing audit')
 source=OUT/'results/II/raw_risk_actions_DEV';closed=OUT/'II/jobs/raw_risk_actions_DEV/CLOSED.json'
 if json.loads(closed.read_text())['status']!='COMPLETE':raise ValueError('closed source required')
 b=Book('II','raw_risk_actions_DEV_float32_audit',max_wall=300,analysis_only=True);status='COMPLETE';error=None
 try:
  dest.mkdir();cfg=json.loads((OUT/'configs/TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
  prior=MaterialPrior(cfg,points);roi=prior.mask;host=np.array(cfg['host_chi']);low=host;high=host+np.array([cfg['inclusion_increment_real_interval'][1],cfg['inclusion_increment_imag_interval'][1]])
  # Half one spacing at each original decimal boundary is sufficient for
  # nearest float32 representation; this is diagnostic, not a changed decoder.
  tolerance=lambda v:.5*np.abs(np.spacing(np.asarray(v,dtype=np.float32))).astype(float)+8*np.finfo(float).eps*np.maximum(1.,abs(v))
  tl=tolerance(low);th=tolerance(high);rows=[];sample=[];spec=json.loads((source/'FROZEN_REQUEST.json').read_text())
  for item in spec['rows']:
   path=Path(item['frozen_mu'])
   if checkpoint_identity(path)!=item['mu_identity']:raise ValueError('frozen action changed')
   with np.load(path,allow_pickle=False) as z:chi=z[item['action']].copy()
   values=np.c_[chi.real,chi.imag];below=np.maximum(low-values,0);above=np.maximum(values-high,0)
   violations=(below>1e-10)|(above>1e-10);substantive=(below>tl)|(above>th)
   outside=np.abs(values[~roi]-host);outside_default=np.any(outside>1e-10,axis=1);outside_sub=np.any(outside>tl,axis=1)
   masks=dict(rectangle_default=np.any(violations[roi],axis=1),rectangle_substantive=np.any(substantive[roi],axis=1),outside_default=outside_default,outside_substantive=outside_sub)
   r=dict(parent=item['parent'],action=item['action'],storage_dtype=str(chi.dtype),default_rectangle_violating_cells=int(masks['rectangle_default'].sum()),
    halfULP_substantive_rectangle_cells=int(masks['rectangle_substantive'].sum()),rounding_only_rectangle_cells=int(np.sum(masks['rectangle_default']&~masks['rectangle_substantive'])),
    max_rectangle_excursion=float(np.maximum(below[roi],above[roi]).max()),max_real_below_host=float(below[roi,0].max()),max_imag_below_host=float(below[roi,1].max()),
    max_real_above_ceiling=float(above[roi,0].max()),max_imag_above_ceiling=float(above[roi,1].max()),
    outside_host_default_cells=int(outside_default.sum()),outside_host_substantive_cells=int(outside_sub.sum()),max_outside_host_component_error=float(outside.max()),
    box_compatible_up_to_float32_rounding=not bool(substantive[roi].any()),ROI_host_compatible_up_to_float32_rounding=not bool(outside_sub.any()),decoder_realizability='NOT_PROVEN',source=str(path))
   rows.append(r);sample.append(dict(**r))
   save(dest/'masks'/item['action']/f"{item['parent']}.npz",roi=roi,real_imag_values=values,below_bound_excursion=below,above_bound_excursion=above,
    lower_halfULP_tolerance=tl,upper_halfULP_tolerance=th,rectangle_default=masks['rectangle_default'],rectangle_substantive=masks['rectangle_substantive'],outside_host_default=outside_default,outside_host_substantive=outside_sub)
  with (dest/'PER_ACTION.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  summaries=[]
  for action in ('prediction_raw','prediction_clipped'):
   selected=[r for r in rows if r['action']==action]
   summaries.append(dict(action=action,parents=len(selected),default_box_violation_parents=sum(r['default_rectangle_violating_cells']>0 for r in selected),
    substantive_box_violation_parents=sum(r['halfULP_substantive_rectangle_cells']>0 for r in selected),substantive_box_violation_cells=sum(r['halfULP_substantive_rectangle_cells'] for r in selected),
    rounding_only_box_cells=sum(r['rounding_only_rectangle_cells'] for r in selected),max_rectangle_excursion=max(r['max_rectangle_excursion'] for r in selected),
    outside_host_substantive_parents=sum(r['outside_host_substantive_cells']>0 for r in selected),max_outside_host_component_error=max(r['max_outside_host_component_error'] for r in selected),decoder_realizability='NOT_PROVEN'))
  result=dict(status='COMPLETE_READONLY_ADDENDUM',source_results=str(source/'RESULTS.json'),source_closed_receipt=str(closed),source_inputs=spec['rows'],summaries=summaries,
   lower_bounds=low,upper_bounds=high,lower_boundary_float32_error=np.array(low,dtype=np.float32).astype(float)-low,upper_boundary_float32_error=np.array(high,dtype=np.float32).astype(float)-high,
   lower_halfULP_tolerance=tl,upper_halfULP_tolerance=th,tolerance_definition='half spacing of nearest float32 boundary plus8eps64*max(1,abs(boundary)); diagnostic only',
   old_default_component_tolerance=1e-10,closed_results_modified=False,materials_changed=False,decoder_realizability='NOT_PROVEN',P1_legal_world='NOT_CLAIMED',new_physics_RHS=0,neural_fits=0)
  write(dest/'RESULTS.json',result)
  lines=['# raw_risk_actions_DEV：只读 float32 边界复核','',
   '原 Full 数组、图像、结果和封账收据保持不变。原检查默认误差 1e-10，远小于 .1/.04 的最近 float32 表示误差；这会把主机常量的表示舍入也列作逐像素违规。冻结材料虽然存为 complex128，常数仍可携带此前 float32 产生的误差。这里按每个原始十进制上下界的半个 float32 spacing，外加很小的 float64 算术余量，单独标记表示舍入可解释与超出该尺度的越界；不据此修改样本。','',
   f'host .1/.04 的 float32 表示误差分别为 {result["lower_boundary_float32_error"].tolist()}；所用下界半 ULP 容差为 {tl.tolist()}。','',
   '|动作|原检查有框越界父对象|半 ULP 后仍有框越界父对象|实质越界格点|仅舍入格点|最大越界|','|---|---:|---:|---:|---:|---:|']
  for r in summaries:lines.append(f"|{r['action']}|{r['default_box_violation_parents']}|{r['substantive_box_violation_parents']}|{r['substantive_box_violation_cells']}|{r['rounding_only_box_cells']}|{r['max_rectangle_excursion']:.9g}|")
  lines += ['', '即使 ROI 外主机和分量框在表示精度下相容，也不证明原始有限维解码器可实现该图像。共享体积分数、共享材料属性、相关场和饱和结构等约束仍未证明，decoder_realizability 始终 NOT_PROVEN，P1_legal_world 仍 NOT_CLAIMED。Full 响应只证明该任意材料动作的数值传播结果，不使其自动成为后验样本。','',
   'PER_ACTION.csv 与 masks/ 保存逐动作误差和分类。新物理、LU、RHS、神经拟合均为 0；成本见 II/jobs/raw_risk_actions_DEV_float32_audit/CLOSED.json。']
  (dest/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n');print(json.dumps(summaries),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

if __name__=='__main__':main()
