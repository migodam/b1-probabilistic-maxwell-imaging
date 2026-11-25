"""Build the self-contained reader page from machine records."""
from collections import defaultdict
from html import escape
import json
from pathlib import Path

from common import HERE


def read(path):
    return json.loads(path.read_text())


def table(heads, rows):
    header="".join("<th scope='col'>"+escape(str(x))+"</th>" for x in heads)
    body="".join("<tr>"+"".join("<td>"+escape(str(x))+"</td>" for x in row)+"</tr>" for row in rows)
    return f"<div class='table'><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>"


def main():
    res=HERE/"results"
    result=read(res/"summary.json")
    names={"actual_vacuum__b0":"实际真空，外部确知真空",
           "actual_vacuum__b0.03125":"实际真空，推断允许弱背景",
           "actual_weak__b0.03125":"实际弱背景，推断允许弱背景"}
    data=[]
    for key,row in result["aggregates"].items():
        side,context=key.split("__",1)
        data.append((side+"²",names[context],row["independent_objects"],
                     f"{row['excluded']:,} / {row['all_charged_centres']:,}",
                     row["witnesses"],row["unresolved"],row["true_centre_retained"],
                     f"{row['neumann_excluded']:,}"))
    main_table=table(("网格","背景合同","对象","严格排除 / 位置检查","物理见证",
                      "未决","真中心保留","方向性 Neumann 排除"),data)
    bro=[]
    for key,row in result["broad"].items():
        name=key.split("__")[1]
        comparison=0 if row["psd_excluded"]==0 else row["psd_excluded"]
        psd=(f"{row['psd_interior_excluded']} / {row['psd_interior_denominator']}（内域）"
             if comparison is None and row.get("psd_interior_denominator") else
             ("仅独立节点复查" if comparison is None else str(comparison)))
        bro.append((name,f"{row['excluded']:,} / {row['denominator']:,}",
                    f"{row['indefinite_best']:,}",psd,
                    ", ".join(map(str,row["state_negative_indices_at_exclusions"]))))
    broad=table(("对象","严格排除","获胜证书 ν > 0","固定权重库 PSD 排除","负指数范围"),bro)
    costs=[]
    for path in sorted(res.glob("sector_*_primary.json")):
        a=read(path)
        spec=path.stem.split("_")
        c=a["costs"]
        costs.append((spec[1]+"²",spec[2],a["clip_classes"],len(a["summary"]),
                      f"{c['setup']['cpu_s']:.3f}",f"{c.get('state_cpu_s',0):.3f}",
                      f"{c.get('online_cpu_s',0):.3f}",f"{a['elapsed']['wall_s']:.2f}"))
    cost_table=table(("网格","形态","剪裁样式","观测数","设定 CPU 秒","状态 CPU 秒",
                      "候选 CPU 秒","总墙钟秒"),costs)
    noise=defaultdict(list)
    for path in sorted((res/"stress").glob("stress_*.json")):
        for key,a in read(path)["rows"].items():
            if a["n_centres"]==1024:
                noise[a["inference_background"],a["noise"]].append(a)
    noise_rows=[]
    for (b,p),rows in sorted(noise.items()):
        noise_rows.append((b,f"{100*p:.0f}%",f"{sum(a['excluded'] for a in rows):,} / {len(rows)*1024:,}",
                           sum(bool(a["witness"]) for a in rows),
                           len(rows)))
    stress=table(("背景上限","噪声","严格排除","物理见证","配对对象"),noise_rows)
    contrasts=[]
    for path in sorted(res.glob("contrast_*.json")):
        for key,row in read(path)["results"].items():
            contrasts.append((key.split("__")[1],str(row["surviving_bin_hull"]),
                              str(row["local_peak_outer_interval"]),
                              f"{row['local_peak_witness']:.3f}" if row["local_peak_witness"] else "未找到",
                              row["global_peak"]["max_excluded"]))
    contrast=table(("对象","存活幅度档包络","逐格独立峰值外界","见证峰值","全域低峰排除上限"),contrasts)
    library=[v["library_all"] for v in result["library"].values() if "library_all" in v]
    incorrect=sum(len(v["wrong_shapes_fully_excluded"]) for v in library)
    by_grid=result["aggregates"]
    small=by_grid["32__actual_weak__b0.03125"]
    large=by_grid["96__actual_weak__b0.03125"]
    primary=f"{small['excluded']:,} / {small['all_charged_centres']:,}"
    native=f"{large['excluded']:,} / {large['all_charged_centres']:,}"
    o=HERE/"index.html"
    html=f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>固定阵列弱背景下的全视场定位</title>
<style>
:root{{--paper:#f7f6ef;--card:#fffdf8;--ink:#1b302c;--a:#17675f;--rule:#d9dfd3;--muted:#4b625c}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.72 system-ui,"PingFang SC",sans-serif}}
main{{max-width:1100px;margin:auto;padding:26px 30px 64px}}h1{{font-size:clamp(2rem,6vw,3.6rem);line-height:1.15;letter-spacing:-.04em;margin:.45em 0}}
h2{{font-size:1.65rem;margin:1.8em 0 .35em}}p{{max-width:81ch}}nav{{display:flex;gap:1rem;flex-wrap:wrap}}a{{color:var(--a);text-underline-offset:.2em}}
summary{{cursor:pointer;color:var(--a);font-weight:650}}figure{{margin:1.1rem 0}}img{{max-width:100%;height:auto;border:1px solid var(--rule);background:white}}
figcaption{{font-size:.9rem;color:var(--muted);max-width:95ch;margin-top:.4rem}}.intro{{border-top:4px solid var(--a);margin-top:.9rem;padding-top:.6rem}}
.table{{overflow-x:auto;background:var(--card);border:1px solid var(--rule);border-radius:.35rem}}table{{width:100%;border-collapse:collapse;font-variant-numeric:tabular-nums;font-size:.88rem;text-align:left}}
th,td{{padding:.5rem .67rem;border-bottom:1px solid var(--rule)}}th{{white-space:nowrap}}tr:last-child td{{border-bottom:0}}
.metrics{{display:grid;gap:.7rem;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));margin:1.1rem 0}}
.metric{{background:var(--card);border:1px solid var(--rule);padding:.65rem .88rem;border-radius:.35rem}}.metric b{{display:block;font-size:1.5rem;color:var(--a)}}
@media(max-width:650px){{main{{padding:15px 15px 32px}}td,th{{padding:.38rem .5rem}}}}
</style></head><body><main>
<nav><a href="../">研究总览</a><a href="README.md">证书与复现</a><a href="results/centres_summary.csv">逐对象 CSV</a><a href="results/summary.json">机器结果</a></nav>
<div class="intro"><h1>固定阵列弱背景下的全视场定位</h1>
<p>400 MHz、2 m × 2 m、7 发 9 收。每个整数格中心都接受检验，边缘候选使用越界剪裁。允许各背景格独立落在 [0, 1/32] 时，两种网格上均获得全视场定位和通过全波复算的物理见证。</p></div>
<div class="metrics"><div class="metric"><b>{primary}</b>32² 实际弱背景，16 件逐中心排除</div>
<div class="metric"><b>{native}</b>96² 实际弱背景，两件逐中心排除</div>
<div class="metric"><b>{small['witness_true_centre']+large['witness_true_centre']} / 18</b>实际弱背景下真中心保留并找到物理见证</div></div>
<h2>四形态库也可定位</h2>
<p>16 件对象共享实盘、环、双夹杂和 Austria 型四种候选形态。检验逐一读入原始复数端口，按位置排除四类的并集。每个候选活动格独立满足轮廓加权的 [3.75, 4.25]。在 {len(library)} 件对象上，{incorrect} 个错误形态的整张图均被排除。</p>
<figure><img src="figures/weak-background-library.png" alt="16 个 32×32 四形态并集空间判定图"><figcaption>绿色表示排除，赭黄表示仍未决，红色表示找到合法见证，空心圆显示评分时读取的真实中心。</figcaption></figure>
<h2>原生 96² 与同赛道对照</h2>
<figure><img src="figures/native-weak-background.png" alt="两件原生对象的材料、反传播和全视场地图"><figcaption>固定噪声比 3%。同一行的材料和见证采用相同颜色范围。第二行对象靠近右边缘，所有越界候选按相同规则裁剪。</figcaption></figure>
{main_table}
<p>实际弱背景下，方向性 Neumann 采用实方向一阶项和全域余项上包，未能通过耦合收缩条件。32² 的 OSM 峰命中 3/16 个中心，前五名含真中心 5/16。OSM 输出排序，三态证书输出对材料类的约束。</p>
<h2>更宽的环允许类使用不定惯性</h2>
<p>核心逐格允许 [3, 5]，自由壳层逐格允许 [0, 5]，背景仍逐格允许 [0, 1/32]。对照固定权重库仅保留正定状态。更宽的连续搜索与任意精度复核在 32² 单独找到至少一个该库遗漏的 PSD 排除。</p>
{broad}
<figure><img src="figures/broad-annulus.png" alt="宽环先验的 32² 与原生 96² 空间结果"></figure>
<p>同一批两件材料的 96²→32² 预测差为原始噪声半径的 5.27 倍和 2.26 倍。比较先对齐端口单位。 <a href="model_transfer/README.md">配对网格转移与阵列核验</a></p>
<h2>改变联合复噪声</h2>
<p>每件对象、每种实际背景固定同一复噪声方向，分别构造 1%、3%、10% 的读数。32² 对四件配对对象全域重做证书。原生用已有距离下界作严格的跨读数转移，结果保存在同一数据目录。</p>
{stress}
<figure><img src="figures/noise-background-phase.png" alt="全视场排除比例随噪声和背景变化"></figure>
<details><summary>位置确定后的幅度和单格峰值</summary>
<p>存活幅度档是同一档位模板的等效系数范围。独立逐格峰值区间允许不同格跨档变化。全域任一单格高值的事件只标记已验算的见证，其余保留未决。</p>{contrast}</details>
<details><summary>实际计算费用</summary><p>按同一形状、两种实际背景的批次计。状态和候选费用是总耗时的组成部分。原生早期批计算中止费用与新运行分账。详细数据见每种形态的 JSON。 <a href="results/development_costs.json">中止和恢复费用</a></p>{cost_table}</details>
<p><a href="README.md">数学约束与复现方法</a> · <a href="freeze.json">冻结合同</a> · <a href="public/acquisition.json">原始观测索引</a> · <a href="results/arb_M32-R01_broad.json">32² 任意精度独立回放</a></p>
</main></body></html>"""
    o.write_text(html)


if __name__=="__main__":
    main()
