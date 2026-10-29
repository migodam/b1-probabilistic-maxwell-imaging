"""Only sealed blind r0 FULL_GN paid endpoint spectra. No physics or truth reads."""
from pathlib import Path
import os,time,json,csv,resource,sys
import numpy as np
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[3]
NAMES=['px','py','rx','ry','contrast','loss','angle']

def main():
 cpu=time.process_time();wall=time.perf_counter();rows=[];vectors=[]
 for case in ['blind_normal_r0','blind_cavity_r0','blind_strong_r0']:
  base=ROOT/'results/online'/case/'FULL_GN';path=base/'endpoint_exact_jacobian.npz';contractpath=base/'ONLINE_INPUT_CONTRACT.json'
  row=dict(case=case,method='FULL_GN',source=str(path.relative_to(ROOT)),status='MISSING')
  if not (base/'ONLINE_SEALED.json').exists():row['reason']='endpoint not sealed';rows.append(row);continue
  if not path.exists() or not contractpath.exists():row['reason']='paid exact endpoint Jacobian or public contract missing';rows.append(row);continue
  with np.load(path) as d:
   need=['z','exact_jacobian','parameter_scales','prior_jacobian']
   missing=[k for k in need if k not in d]
   if missing:row['reason']='missing fields: '+','.join(missing);rows.append(row);continue
   z=d['z'];jc=d['exact_jacobian'];savedscale=d['parameter_scales'];jp=d['prior_jacobian']
  contract=json.loads(contractpath.read_text());spec=contract['spec'];scale=np.asarray(spec['scale']);sigma=float(contract['config']['sigma'])
  if len(z)!=19 or jc.shape[1]!=19 or not np.array_equal(scale,savedscale):row['reason']='unexpected dimension or public/saved scale mismatch';rows.append(row);continue
  J=np.vstack((jc.real,jc.imag))*scale[None,:]/sigma;A=J[:,:12];B=J[:,12:];P=(jp*scale).T@(jp*scale);lamw=P[:12,:12];H=A.T@A+lamw;cross=A.T@B;raw=B.T@B;S=raw-cross.T@np.linalg.solve(H,cross);S=(S+S.T)/2
  svB=np.linalg.svd(B,compute_uv=False);row.update(status='COMPLETE',sigma=sigma,target_names=NAMES,target_scales=scale[12:].tolist(),endpoint_target_coordinates=z[12:].tolist(),B_singular_values_descending=svB.tolist(),B_condition=float(svB[0]/svB[-1]) if svB[-1]>0 else None,wall_prior_precision_scaled=np.diag(lamw).tolist(),target_prior_precision=0,algorithm_damping=0,active_box_coordinates=[NAMES[i] for i in range(7) if z[12+i]<=spec['lower'][12+i]+1e-9 or z[12+i]>=spec['upper'][12+i]-1e-9])
  arrays=dict(z=z,parameter_scales=scale,sigma=np.asarray(sigma),B_scaled_white=B,A_scaled_white=A,target_Gram=raw,target_Schur=S,wall_prior_precision_scaled=lamw,B_singular_values=svB)
  for label,m in [('B_transpose_B',raw),('wall_prior_Schur',S)]:
   vals,V=np.linalg.eigh((m+m.T)/2);vs=V[:,:2].copy()
   for i in range(2):
    if vs[np.argmax(abs(vs[:,i])),i]<0:vs[:,i]*=-1
   row[label]=dict(eigenvalues_ascending=vals.tolist(),matrix_singular_values_descending=np.linalg.svd(m,compute_uv=False).tolist(),eigenvalue_condition=float(vals[-1]/vals[0]) if vals[0]>0 else None,condition_status='POSITIVE' if vals[0]>0 else 'NONPOSITIVE_NUMERICAL_EIGENVALUE',weakest_two_scaled_unit_vectors=vs.T.tolist(),physical_parameter_displacements_per_scaled_unit=(scale[12:,None]*vs).T.tolist(),weakest_two_eigengaps=[float(vals[1]-vals[0]),float(vals[2]-vals[1])])
   arrays[label+'_eigenvalues']=vals;arrays[label+'_eigenvectors']=V;arrays[label+'_weakest_two_scaled']=vs
   for i in range(2):vectors.append(dict(case=case,information=label,weak_order=i+1,eigenvalue=float(vals[i]),**{k:float(vs[j,i]) for j,k in enumerate(NAMES)}))
  np.savez_compressed(OUT/(case+'_spectra.npz'),**arrays);rows.append(row)
 result=dict(scope='Local paid exact Jacobian spectra at sealed blind r0 FULL_GN selected endpoints only; no truth, posterior or global identifiability claim',metric='Real/imag stacked, sigma=.01 white, public parameter scales; B^T B means unprofiled target block, not unwhitened data',vectors='Unit in scaled target coordinates; sign chosen largest absolute component positive; physical displacements separately saved',constraints='Ambient seven-coordinate spectra; active-bound tangent cone is not imposed',cases=rows)
 (OUT/'SPECTRA.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 with (OUT/'WEAK_VECTORS.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['case','information','weak_order','eigenvalue']+NAMES);w.writeheader();w.writerows(vectors)
 text=['# 已封存 r0 成像端点：目标弱方向（离线0PDE）','',result['metric'],'','单位向量按公开目标尺度归一化；符号取最大绝对分量为正。以下为完整七坐标的局部谱，未约化活跃box约束；不是 posterior、置信区间或全局不可辨性。无算法damping，目标prior精度0，墙prior按公开独立Gaussian曲率。','']
 for r in rows:
  text+=['## '+r['case'],'']
  if r['status']!='COMPLETE':text += ['MISSING: '+r.get('reason','unknown'),''];continue
  text += ['来源：`'+r['source']+'`；已检查ONLINE_SEALED。','',f"目标尺度：{r['target_scales']}；活跃box坐标：{r['active_box_coordinates']}。",'', '| 对象 | 升序特征值 | condition |','|---|---|---|']
  for label in ['B_transpose_B','wall_prior_Schur']:
   q=r[label];text += [f"| {label} | "+', '.join(f'{v:.6g}' for v in q['eigenvalues_ascending'])+f" | {q['eigenvalue_condition']:.6g} |"]
  text += ['',f"B奇异值（降序）：{', '.join(f'{v:.6g}' for v in r['B_singular_values_descending'])}；B condition={r['B_condition']:.6g}。",'', '| 对象/弱序 | px | py | rx | ry | contrast | loss | angle |','|---|---:|---:|---:|---:|---:|---:|---:|']
  for label in ['B_transpose_B','wall_prior_Schur']:
   for i,v in enumerate(r[label]['weakest_two_scaled_unit_vectors']):text += [f'| {label}/{i+1} | '+' | '.join(f'{a:+.5f}' for a in v)+' |']
  text += ['']
 text += ['## 解释边界','', '弱向量展示在当前模型、观测布局、尺度和所选端点处哪些shape/contrast/loss组合对数据较弱；墙prior消元允许墙调整后，曲率可能进一步减小。向量中较大的异号分量表示该局部组合可相互抵消，但符号整体可翻转；接近或重合特征值时，应解释弱子空间而非单根向量。这里没有依据真值判断准确性，也没有由小正特征值推断exact fiber。不同端点的谱不是同状态chart误差比较。','', '只有FULL_GN三个已封存r0结果被读取；known-wall分析为可选项，本次不等待、不另求J。']
 (OUT/'REPORT_ZH.md').write_text('\n'.join(text)+'\n')
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
 receipt=dict(physical_cost=dict(full_fe=0,rhs=0),process_cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,peak_rss_bytes=rss,pid=os.getpid(),cases_completed=sum(r['status']=='COMPLETE' for r in rows),read_scope='sealed results/online/*r0/FULL_GN only plus public contract; no inputs/private')
 (OUT/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt));print(json.dumps([{k:r.get(k) for k in ['case','status','B_condition','active_box_coordinates']} for r in rows]))
if __name__=='__main__':main()
