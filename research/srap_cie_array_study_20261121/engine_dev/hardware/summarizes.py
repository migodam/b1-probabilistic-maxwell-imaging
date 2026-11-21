"""Recompute historical matched comparisons and produce the public hardware note."""
from pathlib import Path
import csv,json,statistics as st
from collections import Counter
here=Path(__file__).resolve().parent
data={};gains=[];fields=[]
def table(header,rows):
    return '\n'.join('| '+' | '.join(str(c) for c in row)+' |' for row in [header,['---']*len(header),*rows])
def csvout(name,rows):
    with (here/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for host in ('Mac','XINAN'):
    rows=list(csv.DictReader((here/(host+'_20261120_metrics.csv')).open()))
    cases=sorted({r['case'] for r in rows if r['suite']=='LA_SYNTH' and r['exposed_control'].lower()=='false'})
    assert len(cases)==11 and len(rows)==504
    for view in ('raw','light'):
        rd={(r['case'],r['method']):r for r in rows if r['view']==view and float(r['budget_s'])==30}
        for case in cases:
            vals={k:float(rd[(case,k)]['rel_L2']) for k in ('paper_nie','cie_dual','u11_routed_ap')}
            base=min(vals['cie_dual'],vals['paper_nie'])
            gain=(base-vals['u11_routed_ap'])/base
            gains.append(dict(host=host,case=case,view=view,**vals,reference=base,gain=gain,
                              outcome='W' if gain>0 else 'L' if gain<0 else 'T'))
        for case in sorted({r['case'] for r in rows if r['suite']=='LA_FRESNEL'}):
            vals={k:float(rd[(case,k)]['held_field_err']) for k in ('paper_nie','cie_dual','u11_routed_ap')}
            fields.append(dict(host=host,case=case,view=view,**vals))
for host in ('Mac','XINAN'):
    data[host]={}
    for view in ('raw','light'):
        arr=[r for r in gains if r['host']==host and r['view']==view]
        data[host][view]=[host,view,f"{sum(x['outcome']=='W' for x in arr)}/{sum(x['outcome']=='T' for x in arr)}/{sum(x['outcome']=='L' for x in arr)}",
                          f"{100*st.median(r['gain'] for r in arr):+.3f}%"]
csvout('historical_gain_11.csv',gains)
csvout('historical_fresnel_6.csv',fields)
fres=[]
for case in sorted({r['case'] for r in fields}):
    fr=[case.replace('fresnel_','')]
    for host in ('Mac','XINAN'):
        for arm in ('paper_nie','cie_dual','u11_routed_ap'):
            d={r['view']:r[arm] for r in fields if r['case']==case and r['host']==host}
            fr.append(f"{d['raw']:.3f}/{d['light']:.3f}")
    fres.append(fr)
comp=json.loads((here/'all_old_arrays.json').read_text())
inspection=[]
for folder,key in (('public','Y'),('public','E'),('public','S'),('public','delta'),('public','rx_ids'),('public','valid_mask'),('heldout','Y'),('heldout','E'),('heldout','S'),('heldout','Y_train_opposite_rx')):
    rows=[r for r in comp if (r['folder'],r['key'])==(folder,key)]
    inspection.append([folder,key,len(rows),sum(r['status']=='equal' for r in rows),f"{max(r['relative_F'] for r in rows):.3g}"])
mac=json.loads((here/'mac_fixed_work.json').read_text())
xi=json.loads((here/'xinan_fixed_work.json').read_text())
def med(r,k):
    return st.median(b[k] for b in r['objective_3'])
full=f"""# 两台主机、两份旧记录与相同工作量的审计

旧报告的 **9 胜、2 平、0 负** 来自 Mac 的三阶段、DPCS 秩 15 版本。XINAN 留存的同名旧实验采用四阶段、秩 5，30 秒原始输出为 **2 胜、0 平、9 负**。两份训练包的合成场最大相对差为 4.79e−9，实测训练与留出场逐元素完全相同。当前 XINAN 的九个确认进程各自只有一个操作系统线程。两地差异涉及日程、机器和同时运行的任务数，不能从这两份旧结果估计单独的硬件效应。

## 旧算法的逐组比较

两份原始表均含 504 行。参照原统计代码，剔除已暴露的零旋转 Austria 后，对其余 11 例取 30 秒口径，每一例将 U11 + routed AP 的材料相对 L2 与该例 NIE、Dual-CIE 二者中较低的 L2 做配对。降幅定义为 (基线−U11)/基线，表中中位数对 11 个降幅计算。胜负按原代码的严格大小比较，完全相等计平。

{table(['历史主机','30 s 输出','胜/平/负','L2 中位降幅'],[data[h][v] for h in ('Mac','XINAN') for v in ('raw','light')])}

[Fresnel 和合成的 Mac 旧 CSV](hardware/Mac_20261120_metrics.csv) · [XINAN 旧 CSV](hardware/XINAN_20261120_metrics.csv) · [11 例逐案基线、候选和降幅](hardware/historical_gain_11.csv)

## 旧 Fresnel 的留出场相对误差

每格依次为原始输出 / 轻校正，数值越小越好。两台机器用相同的 Fresnel 阶段参数。实测结果还受单次预算里的物理迭代数、选点和校正的实际工作量影响。

{table(['实测设置','Mac NIE','Mac Dual','Mac U11 AP','XINAN NIE','XINAN Dual','XINAN U11 AP'],fres)}

[六个设置的原精度数值](hardware/historical_fresnel_6.csv)

## 旧启动与运行环境

原会话的运行工具记录保留了实际启动命令。11 月 20 日 XINAN 用 D 盘现行 Python 解释器运行原脚本，合成 16 个并行任务，Fresnel 按脚本降至 8 个。Mac 用项目 Gaussian/.venv_nn 的 Python，合成 6 个并行任务，Fresnel 降至 4 个。两份旧主脚本作逐行比较，只发现两处有实质影响的差异。

| 项目 | Mac 旧脚本 | XINAN 旧脚本 |
| --- | --- | --- |
| 合成傅里叶阶数 | 5,7,10 或 5,8,10 | 3,6,10,12 |
| 合成材料 β 与 DPCS 上限 | 6,2,0.5 或 8,2,0.5；秩 15 | 3,1,0.3,0.1；秩 5 |
| Fresnel 阶段 | 6,10,16；β 6,2,0.5；秩 5 | 同左 |
| Fresnel 最大并行进程数 | 4 | 8 |

两份旧脚本均在导入 NumPy 前将 BLAS/OMP 等线程环境变量设为 1。旧会话没有记录当晚具体 NumPy/SciPy 构建号。11 月 21 日当前解释器路径仍与旧命令对应；XINAN 旧目录和本轮投送目录的 104 个 vendored Python 文件逐字节全部相同。[两条启动记录摘录](hardware/launches.json) · [Mac 旧主脚本](hardware/Mac_20261120_runner.py) · [XINAN 旧主脚本](hardware/XINAN_20261120_runner.py)

## 全量输入对照

逐元素比较 18 个训练 NPZ 和 6 个历史留出 NPZ 的全部数值数组，核对形状、有限值掩码、有限值的最大绝对差与相对 Frobenius 差。所有数组的形状、有限值掩码相同。训练打包中有 6 份接收编号表和一份显式有效性掩码，余下接收协议由完整端口阵列隐式定义。网格坐标、k 与半径均逐位相同。训练场 6 份逐位相同的是全部实测案例，12 份合成训练场呈跨机正演尾差。资料中另有一份可比较的同源留出接收场，其测量值也逐位相同。

{table(['包','数组','数量','逐位相同','最大相对差'],inspection)}

[逐项形状与误差](hardware/all_old_arrays.json) · [直接比较程序](hardware/compare_old_arrays.py)

## 当前主机与等工作量

XINAN 是 Intel i9-14900HX，24 个物理核和 32 个逻辑处理器。只读采样时接通电源，电池为 100%，Windows 使用平衡电源方案。性能计数器给出当时平均约 1853 MHz、处理器时间 61%，另一个时点处理器类的负载为 27%。这两个计数器和时点不同，当次采样无法反推出过去的频率或持续温降。九个二维确认进程均显示 **1 个线程、约 95%–98% 的单核 CPU**，同时另有三维进程工作。使用相同 runner 导入顺序查询运行库，NumPy 和 SciPy 的两份 OpenBLAS 均报告 1 个线程，SciPy FFT 默认 1 个 worker。

Mac mini 使用 Apple M4，10 个物理核。当前 Python {mac['python'].split()[0]}、NumPy {mac['numpy']}、SciPy {mac['scipy']}，两套数学库均链接 Apple Accelerate。线程环境上限为 1，SciPy FFT 报告 1 个 worker；threadpoolctl 在 Mac 上没有枚举出 Accelerate 的动态线程池。XINAN 当前 Python {xi['python'].split()[0]}、NumPy {xi['numpy']}、SciPy {xi['scipy']}，NumPy/SciPy 分别链接 OpenBLAS 0.3.29/0.3.28。

从 180° Austria 旋转 90° 的公开第 16 步取相同材料、入射场、接收算子和 Y，逐字节向两端提供同一压缩输入。在两个当前队列运行期间，各自执行三次 1024×1024 复数 LU、10 个右端和完整物理目标审计，并单独执行 12 对同规格 FFT Green 正/伴随。没有改动电源、进程亲和性或 worker 数量。两台主机目标值之差为 {abs(mac['objective_3'][0]['phi']-xi['objective_3'][0]['phi']):.2g}。这次短测给出局部吞吐，不能把全部反演按同一倍数换算。

{table(['工作量','Mac 墙钟','XINAN 墙钟','XINAN/Mac'],[
['物理审计，3 次中位数',f"{med(mac,'seconds'):.4f} s",f"{med(xi,'seconds'):.4f} s",f"{med(xi,'seconds')/med(mac,'seconds'):.2f}×"],
['12 对 Green 正/伴随',f"{mac['green_pairs_12']['seconds']:.4f} s",f"{xi['green_pairs_12']['seconds']:.4f} s",f"{xi['green_pairs_12']['seconds']/mac['green_pairs_12']['seconds']:.2f}×"]])}

当前按主机分配的独立确认比较，应以相同案例、相同布局、同机同预算的双臂为准；历史的 9/11 只对应其原 Mac 日程及旧布局。[Mac 实测](hardware/mac_fixed_work.json) · [XINAN 实测](hardware/xinan_fixed_work.json) · [相同工作量代码](hardware/fixed_work.py)
"""
(here.parent/'HARDWARE.md').write_text(full)
