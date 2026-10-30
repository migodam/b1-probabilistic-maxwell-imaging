"""Build the existing-repository Pages entry directly from closed receipts."""
from html import escape
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPOSITORY = "https://github.com/migodam/b1-probabilistic-maxwell-imaging"
BASE = "research/research_z_alg_cuda_pipeline/"


def main():
    data = json.loads((ROOT / "RESULTS_SUMMARY.json").read_text())
    runs = data["endtoend_runs"]
    table = "".join(f'<tr><td>{r["lanes"]} lane'+(' + hybrid SVD' if 'hybrid' in r['run'] else ' / CPU SVD')+
                    f'</td><td>{r["wall_s"]:.3f} s</td><td>{r["optimizer_evaluations"]}</td><td>{r["speedup_vs_one_lane"]:.3f}×</td></tr>' for r in runs)
    links = "".join(f'<a href="{REPOSITORY}/blob/main/{BASE}{path}">{label}</a>' for path, label in [
        ("RESEARCH_REPORT_ZH.md", "中文完整报告"), ("OPTIMIZED_OPERATOR_CATALOG.md", "可复用算子目录"),
        ("REPRODUCE.md", "接入与复现"), ("RESULTS_SUMMARY.json", "实际结果 JSON"),
        ("COST_LEDGER.csv", "费用账本"), ("src/device.py", "CUDA 常驻算子"),
        ("src/executor.py", "任务 wrapper"), ("src/hybrid_svd.py", "混合 SVD"), ("FAILURES.json", "失败记录")])
    html = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Maxwell CUDA 执行流水线与现成算子复用</title><style>
:root{{--ink:#142a39;--muted:#506775;--line:#d8e3e9;--accent:#16674e}}*{{box-sizing:border-box}}body{{margin:0;background:#f3f7f8;color:var(--ink);font:16px/1.75 system-ui,-apple-system,sans-serif}}main{{max-width:1160px;margin:auto;padding:42px 24px 64px}}header{{border-bottom:2px solid var(--line);padding-bottom:26px}}.tag{{color:var(--accent);font-size:14px;letter-spacing:.08em}}h1{{font-size:clamp(27px,4vw,42px);line-height:1.25;margin:12px 0 18px}}h2{{font-size:25px;margin-top:42px}}h3{{font-size:19px}}p{{max-width:940px}}.muted{{color:var(--muted)}}.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:28px 0}}.metric,section.box{{background:white;border:1px solid var(--line);border-radius:12px;padding:18px}}.metric b{{display:block;font-size:29px;color:var(--accent)}}.metric span{{font-size:13px;color:var(--muted)}}img{{width:100%;height:auto;background:white;border:1px solid var(--line);border-radius:10px}}figure{{margin:24px 0}}figcaption{{font-size:13px;color:var(--muted);margin-top:8px}}table{{border-collapse:collapse;background:white;width:100%;font-variant-numeric:tabular-nums}}td,th{{padding:12px 15px;text-align:left;border-bottom:1px solid var(--line)}}th{{background:#e9f0f3}}a{{color:#14618b}}nav{{display:flex;flex-wrap:wrap;gap:10px;margin:18px 0}}nav a{{padding:7px 12px;border:1px solid var(--line);background:white;border-radius:6px;text-decoration:none;font-size:14px}}code{{background:#e6eef1;padding:2px 5px;border-radius:4px;font-size:.9em}}pre{{overflow:auto;background:#162e3b;color:#eaf2f4;border-radius:8px;padding:20px;font-size:14px;line-height:1.6}}.columns{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}footer{{margin-top:42px;padding-top:22px;border-top:1px solid var(--line);font-size:13px;color:var(--muted)}}@media(max-width:700px){{.metrics,.columns{{grid-template-columns:1fr 1fr}}main{{padding:25px 14px}}td,th{{padding:9px 6px;font-size:13px}}}}
</style></head><body><main><header><div class="tag">MAXWELL · CUDA EXECUTION · OPERATOR REUSE</div><h1>把可并行的计算持续交给 GPU</h1>
<p>在完整 576 单元非线性 Maxwell 反演中复用成熟算子：常驻物理、私有 stream、连续任务队列，以及不阻塞 CPU 的混合大 SVD。</p>
<p class="muted">RTX 4060 Laptop / i9-14900HX · complex128 / float64 · 已公开观测的计算工程实验 · 无新 GT recovery 主张</p></header>
<div class="metrics"><div class="metric"><b>43.0%</b><span>同一批次总 wall 下降</span></div><div class="metric"><b>56.715 s</b><span>12 个完整多尺度任务</span></div><div class="metric"><b>631</b><span>各配置相同目标评估数</span></div><div class="metric"><b>635 LU</b><span>各配置相同物理求解量</span></div></div>
<nav>{links}<a href="{REPOSITORY}/blob/main/AGENTS.md">仓库 AGENTS.md</a></nav>
<section><h2>先看完整任务，而不是只看小内核</h2><p>三个匿名公开目标，各运行 Native TRF、Reciprocal TRF、data L-BFGS-B 和 GlobalDual+data correction。共同 8→12→24 材料参数化，物理始终为 576 阶完整 VIE。所有配置均正常返回 12 个任务，保持 631 目标评估、635 LU、26,880 complex RHS。正常返回不等于全部恢复材料。</p>
<table><thead><tr><th>执行路径</th><th>总 wall</th><th>目标评估</th><th>相对单 lane 吞吐</th></tr></thead><tbody>{table}</tbody></table>
<figure><img src="figures/end_to_end_queue.png" alt="完整任务费用和四 lane 连续任务时间线"><figcaption>下方是 CPU+CUDA 的任务区间，并非 GPU kernel 或 SM occupancy trace。总 wall 含共享 GPU setup、任务和写盘；不含 Python/Torch 导入及 SSH 传输。</figcaption></figure></section>
<section><h2>怎样减少人为等待</h2><div class="columns"><div class="box"><h3>物理常驻与一次返回</h3><p>G/E/S 只上传一次。每个材料的全部源共享 LU，forward、adjoint、Gram、正则和粗网格乘法在 GPU 连续完成。值与梯度打包到 pinned buffer，只等待自己的 completion event。</p></div><div class="box"><h3>有空位就用，忙时继续 CPU</h3><p>独立任务各用私有 stream。大 SVD 最多一个 GPU slot；slot 已被占用时，其他任务立即调用原 CPU SVD。64/144 维 SVD 留在 CPU，576 维允许使用现成 cuSOLVER。</p></div></div>
<p>混合执行中 318 次 SVD 有 46 次在 GPU。并发期间 GPU 大 SVD 单次平均并没有更快，但分担 CPU 通道、同时推进 ready work，使整批任务缩短。并行子任务 wall 不能相加代替整批 wall。</p></section>
<section><h2>哪些工作适合批量交给 CUDA</h2><figure><img src="figures/cuda_throughput.png" alt="目标梯度和候选评分真实吞吐"><figcaption>数据梯度与 GlobalDual 的 4-stream 吞吐分别约 130.6 / 123.5 次每秒；batch=16 全非线性候选评分约 221.9 个材料每秒。候选评分仅含已驻留输入的正演，不能解释成完整反演速度。</figcaption></figure>
<p>多初值、多目标可作为独立任务；同一材料的多 illumination 应作为一个 LU 的多 RHS。普通 stream 不对应固定 SM 分区，本轮没有改驱动、安装 MPS 或新调度服务。</p></section>
<section><h2>已有的好算子，先复用再开发</h2><p><a href="{REPOSITORY}/blob/main/{BASE}OPTIMIZED_OPERATOR_CATALOG.md">算子目录</a>列出 A9 的 DenseDDA CUDA LU、多 RHS；A17 压缩 endpoint 的批量 Cholesky/Woodbury；Route A/B 的 FFT 传播；R2/R9/R10 的固定块/宿主/形状复用，并区分实际采用与仅审计。向量 3D、固定 ROI 和 scalar-TM 的物理条件不能直接互换。</p>
<figure><img src="figures/svd_operator_selection.png" alt="CPU和现成cuSOLVER SVD算子尺寸选择"><figcaption>实际公开 Jacobian 的尺寸选型实验，计入必要 GPU 上传/下载；单次且含冷 API 费用，不作跨机器统计保证。</figcaption></figure>
<p>本机项目根和本仓库的 <code>AGENTS.md</code> 已加入复用指引。其他 session 可以直接采用相同任务 wrapper；没有向其他聊天发消息或改动它们的任务。</p></section>
<section class="box"><h2>可直接运行</h2><pre>D:\\python\\python.exe -m src.executor run manifest_4_hybrid.json</pre><p>输入是公开 G/E0/S/Y；推荐入口使用 4 lane + 非阻塞混合 SVD。具体路径、四种方法、硬件、计时窗与失败记录见<a href="{REPOSITORY}/blob/main/{BASE}REPRODUCE.md">复现说明</a>。</p></section>
<section><h2>真实失败和结论边界</h2><p>默认库 CUDA Graph 捕获失败；换成 cuSOLVER 后成功，却没有成为最快的批量评分路径。强制 cuSOLVER 也使 batch=16 变慢，因此默认保留库的原选择策略。现有 profiler 仅支持 CPU，SM occupancy 与 kernel overlap 为 NOT_MEASURED。</p>
<p>本成果是成熟物理/线性代数原语的数据驻留与资源调度组合。没有新 Maxwell 定理、概率解集、严格分支覆盖、全局收敛或新 GT recovery 结论；旧 R2/R3/R4 研究保持冻结。</p></section>
<footer>2026-10-30 · <a href="../research_z_alg_r4_diversity/">R4-DIVERSITY</a> · <a href="https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/{BASE}RESEARCH_REPORT_ZH.md">匿名 raw Markdown</a> · 所有图、数字与账本来自实际完成回执。</footer></main></body></html>'''
    (ROOT / "index.html").write_text(html, encoding="utf-8")
    print(ROOT / "index.html")


if __name__ == "__main__":
    main()
