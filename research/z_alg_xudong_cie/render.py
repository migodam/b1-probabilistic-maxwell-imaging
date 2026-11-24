"""Reader-facing report using the existing Maxwell project stylesheet."""
import csv
from html import escape
import json
import os
from pathlib import Path
import re
import markdown

BASE=Path(__file__).resolve().parent
DIR="research/z_alg_xudong_cie"
REPO="https://github.com/migodam/b1-probabilistic-maxwell-imaging"
RELEASE=REPO+"/releases/tag/z-alg-xudong-cie-20261124"
ZIP=REPO+"/releases/download/z-alg-xudong-cie-20261124/Z_ALG_Xudong_CIE_Research.zip"
PAGES={
 "RESEARCH_REPORT_ZH.md":("index.html","研究结论"),
 "NEXT_PRO_HANDOFF.md":("handoff.html","ChatGPT 交接"),
 "CURRENT_THEORY.md":("current-theory.html","电流与材料"),
 "MINIMUM_PEAK_THEORY.md":("peak-theory.html","峰值与面积"),
 "pruning/THEOREM.md":("pruning/support-theory.html","支撑排除证明"),
 "pruning/OUTSIDE_CONTRAST.md":("pruning/outside.html","区外材料"),
 "diagnosis/REPORT_ZH.md":("diagnosis/index.html","四个开发病例"),
 "minpeak_native/REPORT_ZH.md":("minpeak_native/index.html","六件原生材料"),
 "minpeak_maxstruct/REPORT_ZH.md":("minpeak_maxstruct/index.html","七件存档材料"),
 "area/README.md":("area/index.html","面积下界"),
 "current_physics/CANCELLATION_ZH.md":("current_physics/index.html","守恒相消"),
 "chambers/README.md":("chambers/index.html","暗室数据"),
 "cross_cie/REPORT_ZH.md":("cross_cie/index.html","守恒更新与仅重启的对照"),
 "pruning/REPORT_ZH.md":("pruning/index.html","候选材料位置"),
 "pruning/site_control/RESULTS_ZH.md":("pruning/site_control/index.html","四个缺陷的位置盲测"),
 "pruning/CONSTITUTIVE_MARGIN.md":("pruning/constitutive.html","辅助电流的最小本构误差"),
 "pruning/RECIPROCITY.md":("pruning/reciprocity.html","能量与互易性反例"),
 "mechanism/REPORT_ZH.md":("mechanism/index.html","逐像素光学相位的独立试验"),
 "holdout/projection/README.md":("holdout/projection/index.html","粗网格能承载多少材料细节"),
 "holdout/README.md":("holdout/index.html","新材料和新端口"),
 "NEON_GUIDANCE_ZH.md":("guidance.html","外部研究意见与检验"),
 "PHYSICS_AND_SOURCES.md":("sources.html","物理与原始文献"),
 "ARCHIVE_REPLAY.md":("archive.html","研究包和原始实测"),
 "REPRODUCE.md":("reproduce.html","复算说明")}
NAMES={"dual":"原 Dual-CIE","reset":"仅重启共轭梯度","diag":"逐源光学修正","full":"全矩阵光学与扇区","pixel":"逐像素光学修正",
       "multistart":"三个普通初值","matched_boxes":"分支全算","pruned":"证书剪枝（单件）",
       "pruned_reuse":"证书剪枝（复用几何）","dual_then_restart":"Dual＋普通第二初值","dual_then_pixel":"Dual＋逐像素第二路线"}
def relative(p,output):
    return os.path.relpath(BASE/p,(BASE/output).parent)
def template(output,title,body):
    def l(path):return escape(relative(path,output))
    nav=[("index.html","中文报告"),("results.html","逐件结果"),("reproduce.html","复算"),("handoff.html","ChatGPT 交接")]
    tabs=" ".join(f'<a href="{l(link)}"'+(' aria-current="page"' if output==link else '')+f'>{label}</a>' for link,label in nav)
    return f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)} · Maxwell 逆散射研究</title><link rel="stylesheet" href="{l("../z_alg_optical_inertia/research.css")}">
<style>details{{margin:1.1rem 0}}summary{{cursor:pointer;color:var(--accent);font-weight:650}}.math-display{{overflow-x:auto}}figure img,article img{{max-width:100%;height:auto}}details table{{font-variant-numeric:tabular-nums}}.table-scroll{{overflow:auto}}</style>
<script>window.MathJax={{tex:{{inlineMath:[['$','$'],['\\\\(','\\\\)']]}}}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js"></script></head><body>
<a class="skip-link" href="#content">跳到正文</a><header class="site-header"><div class="site-shell"><div class="project-links">
<a href="{l("../../")}">B1 总项目</a><a href="{l("../experiment_handbook/")}">实验手册</a>
<a href="{l("../z_alg_optical_inertia/ip_submission_v2/manuscript_ip.pdf")}">Paper 1 英文稿</a><a class="download" href="{ZIP}">下载本轮研究包</a></div><nav>{tabs}</nav></div></header>
<main class="site-shell" id="content"><article>{body}</article><footer><a href="{REPO}/tree/main/{DIR}">源码与数据说明</a> · <a href="{RELEASE}">版本与附件</a></footer></main></body></html>'''
def render_md(src,dst):
    text=(BASE/src).read_text()
    saved=[]
    def protect(match):
        idx=len(saved);saved.append(match.group(0));return f'ZALGFORMULATOKEN{idx}Z'
    text=re.sub(r'\$\$.*?\$\$|\\\[.*?\\\]|\\\(.*?\\\)|(?<!\$)\$(?!\$)[^\n$]+\$',protect,text,flags=re.S)
    content=markdown.markdown(text,extensions=["tables","fenced_code","sane_lists","toc"])
    for j,expr in enumerate(saved):
        content=content.replace(f'ZALGFORMULATOKEN{j}Z',escape(expr))
    def link(m):
        url=m.group(1)
        if "://" in url or url.startswith(("#","/","mailto:")):return m.group(0)
        path,sep,anchor=url.partition("#")
        absolute=(BASE/src).parent/path
        if not path:return m.group(0)
        try:
            where=absolute.resolve().relative_to(BASE)
        except ValueError:
            if path.endswith((".py",".md",".json")):
                resolved=absolute.resolve()
                project=BASE.parents[1]
                if resolved.is_relative_to(project):
                    if "z_alg_optical_inertia" in resolved.parts:
                        url=REPO+"/releases/tag/z-alg-optical-inertia-validation-v1"
                    elif "srap_cie_array_study_20261121" in resolved.parts:
                        url=REPO+"/releases/tag/srap-cie-array-study-20261121"
                    else:
                        url=RELEASE
                    return f'href="{escape(url,quote=True)}"'
            return m.group(0)
        if str(where) in PAGES:
            url=relative(PAGES[str(where)][0],dst)+sep+anchor
        elif where.suffix in (".md",".py",".csv",".json"):
            on_site=(absolute.is_file() and absolute.stat().st_size<=384*1024
                     and (len(where.parts)<=2 or where.parts[:2] in (("pruning","site_control"),("holdout","projection"))))
            url=f"{REPO}/blob/main/{DIR}/{where}"+sep+anchor if on_site else RELEASE
        elif where.suffix in (".npz",".log") or absolute.is_dir():
            url=RELEASE
        return f'href="{escape(url,quote=True)}"'
    return re.sub(r'href="([^"]*)"',link,content).replace("<table>",'<div class="table-scroll"><table>').replace("</table>","</table></div>")
def rows_table(headers,rows):
    out='<div class="table-scroll"><table><thead><tr>'+''.join(f'<th>{escape(x)}</th>' for x in headers)+'</tr></thead><tbody>'
    for row in rows:
        out+="<tr>"+"".join(f"<td>{escape(str(v))}</td>" for v in row)+"</tr>"
    return out+"</tbody></table></div>"
def fmt(v):
    return f"{float(v):.3f}" if v not in ("",None) else "未决"
def per_object():
    rows=list(csv.DictReader((BASE/"benchmark/per_observation.csv").open()))
    lab=list(csv.DictReader((BASE/"benchmark/laboratory_all.csv").open()))
    groups=sorted({r["object"] for r in rows})
    parts=['<h1>同一物体、同一观测，各方法交出的材料</h1><p>材料误差由原生 96² 的参考图计算。目标内、空气背景的两项均以真材料全图的范数归一，平方和等于全图相对误差的平方。IoU 越大越好，其他误差越小越好。未见端口使用另一套收、发位置。墙钟包含训练端口的全波选图和证书计算。</p>']
    parts.append('<h2>三件暗室实物</h2><p>未见发射的复场误差按真实有效通道计。全图、内部、空气和 IoU 列以实验室给出的标称几何和材料值作参考，空气列为绝对均方根误差。</p>')
    for case in sorted({r['observation'] for r in lab}):
        selected=[r for r in lab if r["observation"]==case]
        parts.append(f'<details><summary>{escape(case)}</summary>'+
             rows_table(["方法","训练复场","未见发射","标称全图","标称内部","空气 RMS","标称 IoU","墙钟（秒）"],
                        [[NAMES.get(r["arm"],r["arm"]),*[fmt(r.get(k)) for k in
                         ("train","held","nominal_full_grid_relative_l2","nominal_interior_relative_l2",
                          "nominal_air_RMSE","nominal_support_IoU","wall")]] for r in selected])+'</details>')
    parts.append('<h2>十二件独立合成材料</h2>')
    for oid in groups:
        parts.append(f'<details id="{oid}"><summary>{oid}</summary>')
        for ap in ("pw180_10x20","pw360_10x20"):
            for noise in ("n02","n10"):
                selected=[r for r in rows if (r["object"],r["aperture"],r["noise"])==(oid,ap,noise)]
                idx={k:i for i,k in enumerate(NAMES)}
                selected.sort(key=lambda r:idx.get(r["arm"],100))
                parts.append(f'<h3>{ap[2:5]}° 孔径 · {int(noise[1:])}% 噪声</h3>')
                parts.append(rows_table(["方法","全图","目标内","空气","半峰 IoU","训练复场","未见复场","墙钟（秒）"],
                    [[NAMES.get(r["arm"],r["arm"]),*[fmt(r[k]) for k in ("L2","inside","air","IoU","train","held","wall")]] for r in selected]))
        parts.append("</details>")
    return '\n'.join(parts)
if __name__=="__main__":
    for src,(dest,title) in PAGES.items():
        (BASE/dest).write_text(template(dest,title,render_md(src,dest)))
    (BASE/"results.html").write_text(template("results.html","逐件结果",per_object()))
