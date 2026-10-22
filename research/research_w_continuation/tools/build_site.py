"""Static public reading page from completed cached evidence, no physics."""
from pathlib import Path
import json, html
root=Path(__file__).resolve().parents[1]
gate=json.loads((root/'GATE_DECISION.json').read_text())
q=json.loads((root/'P0_Q_QUALIFICATION.json').read_text())
cost=json.loads((root/'COST_SUMMARY.json').read_text())['stages']
repo='https://github.com/migodam/b1-probabilistic-maxwell-imaging'
base=repo+'/blob/main/research/research_w_continuation/'
def doc(path,label):return f'<a href="{base}{path}">{html.escape(label)}</a>'
def figure(name,title,caption):
    return f'<figure><a href="figures/{name}.png"><img loading="lazy" src="figures/{name}.png" alt="{html.escape(title)}"></a><figcaption><strong>{html.escape(title)}</strong> — {html.escape(caption)}</figcaption></figure>'
table=''.join(f'<tr><td>{s}</td><td>{r["reserved_complex_rhs"]}</td><td>{r["reserved_factorizations"]}</td><td>{r["charged_cpu_seconds"]:.2f}</td><td>{r["measured_job_wall_seconds"]:.2f}</td></tr>' for s,r in cost.items())
document=f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Research W continuation｜物理算准了，传播歧义尚未得到见证</title>
<meta name="description" content="高阶二维Maxwell正演资格通过；三组四方法搜索未发现合格传播歧义见证。真实数组、负结果、费用及复现。">
<style>:root{{--ink:#17283b;--muted:#58677a;--line:#dce3eb;--blue:#176b9c;--bg:#f5f7fa}}*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:17px/1.85 system-ui,-apple-system,"PingFang SC",sans-serif}}main{{max-width:1080px;margin:auto;padding:45px 24px 80px}}h1{{font-size:clamp(29px,4.5vw,47px);line-height:1.35;letter-spacing:-.02em;max-width:950px}}h2{{margin:42px 0 12px;font-size:25px}}h3{{font-size:19px}}p{{max-width:950px}}a{{color:var(--blue);text-underline-offset:4px}}.eyebrow{{font-size:14px;color:var(--muted);letter-spacing:.08em}}.lead{{font-size:20px;max-width:930px}}.status{{display:flex;gap:10px;flex-wrap:wrap}}.chip{{padding:5px 11px;border-radius:5px;background:#e0f0e8;color:#176443;font-size:14px}}.chip.hold{{background:#fff0d4;color:#80551d}}nav{{display:flex;gap:9px 22px;flex-wrap:wrap;margin:24px 0}}section.card,.note{{background:white;border:1px solid var(--line);border-radius:10px;padding:22px 25px;margin:22px 0}}.note{{border-left:4px solid #bd7824}}figure{{margin:28px 0;background:white;border:1px solid var(--line);border-radius:10px;overflow:hidden}}img{{display:block;width:100%;height:auto}}figcaption{{padding:14px 19px;color:var(--muted);font-size:14px}}table{{border-collapse:collapse;width:100%;font-size:15px}}th,td{{text-align:left;border-bottom:1px solid var(--line);padding:10px}}.table{{overflow-x:auto}}small,footer{{color:var(--muted);font-size:13px}}li{{margin:9px 0}}summary{{cursor:pointer;color:var(--blue)}}code{{font-size:.85em;word-break:break-word}}footer{{border-top:1px solid var(--line);padding-top:20px;margin-top:45px}}</style></head><body><main>
<div class="eyebrow">RESEARCH W · CONTINUATION · 2026-10-23 · 二维 TM / CPU</div>
<h1>物理算准了，<br>传播歧义尚未得到见证</h1>
<p class="lead">先确认两种材料响应的差别不是网格误差，再问未知墙体是否会留下影响目标成像的传播多解。本轮完成高精度前向修复与有限搜索，没有启动概率或递进成像。</p>
<div class="status"><span class="chip">P0-Q：离散全波资格通过</span><span class="chip hold">P0-B：分配的搜索额度内未找到</span><span class="chip hold">joint / P1：未运行</span></div>
<nav>{doc('RESEARCH_W_CONTINUATION_REPORT_ZH.md','完整中文报告')}{doc('P0_B_CALIBRATION_RESULTS.csv','逐候选数据')}{doc('REPRODUCE.md','代码与复现')}<a href="{repo}/releases/tag/research-w-continuation-20261023">原始复场数组 Release</a><a href="https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_w_continuation/RESEARCH_W_CONTINUATION_REPORT_ZH.md">匿名 raw 报告</a></nav>
<section class="card"><h2 style="margin-top:0">三个直接结论</h2><ol><li><strong>旧误差主要是空间离散问题。</strong>只提升Green积分阶数几乎无效；贴合材料界面并采用高阶场表示后，独立线源圆柱相对误差降至0.0266%。没有拟合自由复增益。</li><li><strong>原始墙体的完整观测达到既定精度。</strong>两seed、两带、cal及固定目标共8项最细层级差为0.00653–0.03322，小于原0.1；PML、圆柱、伴随和互易另有检验。这是有限的离散资格，不是连续PDE认证。</li><li><strong>尚未进入“数据近等价”的物理区域。</strong>3个context×4方法中最好cal距离为<strong>{gate['minimum_observed_dcal']:.4f}</strong>，门槛为1。没有任何合格歧义对，也不能从未找到推断唯一性。</li></ol></section>
<h2>先验证数值模型，而非只验证矩阵恒等式</h2><p>升级VIE与独立高阶FEM使用相同连续材料、源、接收器及噪声单位。双单元平均VIE的负结果也保留。最终选用界面相容P3 FEM、解析入射和PML；保持全部二维多次散射。</p>
{figure('01_cylinder_accuracy','独立Mie参考与受控精度改动','图中的所有误差来自真实数值数组；旧VIE结果与新FEM控制分别标明。')}
{figure('02_wall_convergence','原两墙体的八项完整观测收敛','完整白化范数，不除以通道数平方根；相邻层级差不是严格连续误差上界。')}
{figure('11_saved_task_total_field','实际保存的总场幅值与相位','原两seed，中心频率、第一源；积分点确定性抽样，不生成任何新物理场。')}
<details><summary>查看核积分、PML和参考复场</summary>{figure('03_kernel_quadrature','近场与自单元积分','原生产self使用解析积分；虚线self GL只是诊断，不是原求解器实际误差。')}{figure('05_pml_extraction','PML与接收提取检查','外域变化包含重网格效应；不同敏感度不冒充严格可加误差预算。')}{figure('04_known_world_and_field','参考墙体与实际复场','这是已知世界的正演结果，不是逆成像成功图。')}</details>
<h2>四方法公平搜索：所有候选都离相容门槛很远</h2><p>Full-J、同材料尺度随机、双初始化FWI和Task-aware使用同一全域材料坐标、约束、先验与正确器权限。提案和接受不读取truth。这里是短预算局部搜索，并非充分收敛的多启动FWI。</p>
{figure('06_p0b_calibration_and_rhs','cal距离与独立运行RHS归因','共享setup实际支付一次；归因柱不能相加当作总费用，也不是端到端加速结果。')}
{figure('07_p0b_offline_material_comparison','已知参考、共同BP初值与在线选定端点','truth只用于输出保存后的评分；图中较差拟合明确保留，未包装成恢复成功。')}
{figure('09_p0b_calibration_cost_scatter','误差—RHS费用散点','RHS为独立运行归因，非完整wall-time；所有点仍未达到cal距离1。')}
<div class="note"><strong>为何只有3/4个context？</strong>首次搜索发现正确器回退错误，中止并保留760 RHS费用。修正版完成三组后剩480搜索RHS，不足完整800 RHS四方法预留；没有挪用端点/验证额度。全局预算尚有余额，不等于可重新扩展冻结搜索合同。</div>
<div class="note"><strong>另有明确的表示实施缺口。</strong>原0.1λ最小特征没有被48×32像素及投影严格执行；box、passivity和面积合法不等于所有几何先验合法。当前候选不能用于正向物理见证，也不能据搜索失败排除原合法材料类中的多解。未事后放宽原限制。</div>
<h2>目标响应与joint必须分开</h2><p>所有候选只计算了原24任务中的固定目标8。原未贴合圆形目标的网格细化差0.198–0.270，不合格并保留。三个在线选定端点补做圆周贴合后，256→512段差为0.000644–0.000923，求解残差通过；这个数值修复不能把远大于1的cal距离变成传播歧义。</p><p>相同目标的质心差为0。未知目标joint补偿、全部24任务、最终见证扰动、posterior、A17新ROM、NN及三维均未运行。</p>
{figure('08_p0b_fixed_probe8','原未贴合目标网格输出','仅作原始诊断，三个端点的贴合复核另列；没有合格物理见证。')}
{figure('10_p0b_task_refinement_comparison','同一端点的界面处理复核','贴合前与贴合后比较口径不同，均为经验一致性；保留极瘦单元诊断及全部安全检查，不是连续PDE误差界。')}
<h2>费用与失败</h2><div class="table"><table><thead><tr><th>新阶段</th><th>预约复RHS</th><th>预约分解</th><th>计费CPU秒</th><th>作业wall秒</th></tr></thead><tbody>{table}</tbody></table></div><p><small>这是报告快照；发布尾项收据单列。包含失败和重跑，预约与完成数在机器账本分别列出。Q计费含60秒未仪表化准备分配，不能称为精确实测CPU。旧P0费用不混入。PML审核峰值2.160GB超过原1.5GB工程目标，未发生系统OOM；修正版搜索峰值约1.059GB。</small></p>
<nav>{doc('COST_SUMMARY.json','完整费用')}{doc('PHYSICS_LEDGER.jsonl','物理账本')}{doc('FAILURE_LEDGER.jsonl','失败记录')}{doc('P0_B_PHYSICS_VALIDATION.md','候选数值审核')}</nav>
<h2>下一轮的科学问题</h2><p>应先让普通受约束FWI在原观测和声明先验下达到相容尺度，区分优化失败、表示限制和真实传播多解。只有稳健的joint见证才支持低维概率传播试验；“递进”还需在直接多假设方法之上取得独立收益。本轮没有证据授予这些算法GO。</p>
<nav>{doc('SCIENTIFIC_DECISION.md','科学决策')}{doc('P0_C_IDENTIFIABILITY_PLAN.md','可辨识性方案')}{doc('P1_A_POSTERIOR_PLAN.md','未来概率参考')}{doc('P1_B_PROGRESSIVE_PLAN.md','未来2×2消融')}</nav>
<footer>原P0保留于commit db7146cae2d74ba231577d42f225dbfc9a389268。本页只报告实际完成的二维离散全波实验；未来方案不写作成功结果。源码、配置、种子、复观测、复场、费用和失败均公开于同一仓库。</footer>
</main></body></html>'''
(root/'index.html').write_text(document)
print(json.dumps(dict(page='index.html',physics_actions=0,scientific_status=gate['P0_B_calibration'])))
