"""Export the own R7 report and evidence archive. No scientific execution."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import html
import json
import re
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
RAW = 'https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r7/'
RELEASE = 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r7-20261104'


def markdown_report(source):
    # Retain the project's minimal renderer for prose; handle formula fences
    # separately so line breaks and literal operators survive HTML export.
    import sys
    sys.path.insert(0, str(ROOT.parent/'research_max_struct_r1/experiments'))
    from export_delivery import markdown_subset
    pieces = re.split(r'^```[^\n]*\n(.*?)^```\s*$', source, flags=re.M | re.S)
    body = '\n'.join('<pre><code>'+html.escape(part.rstrip())+'</code></pre>'
                     if i % 2 else markdown_subset(part)
                     for i, part in enumerate(pieces))
    body = re.sub(r'<p>(<figure>.*?</figure>)</p>', r'\1', body)
    return re.sub(r'href="([^":]+\.md)"', lambda m: 'href="'+RAW+m[1]+'"', body)


def main():
    tick,cpu=time.perf_counter(),time.process_time()
    body=markdown_report((ROOT/'RESEARCH_REPORT_ZH.md').read_text())
    css="body{margin:0;background:#f3f5f7;color:#172230;font:17px/1.8 system-ui,-apple-system,sans-serif}main{max-width:1200px;margin:32px auto;background:white;padding:40px;border-radius:12px}h1{font-size:30px;line-height:1.4}h2{margin-top:40px;border-top:1px solid #dde3e8;padding-top:18px}a{color:#075b93}nav{display:flex;gap:16px;flex-wrap:wrap}table{border-collapse:collapse;width:100%;font-size:14px;display:block;overflow:auto}th,td{border:1px solid #d9e1e8;padding:9px;text-align:left}th{background:#eaf1f6}img{max-width:100%;height:auto}figure{margin:24px 0}code{background:#eef2f4;overflow-wrap:anywhere}pre{padding:18px;overflow:auto;background:#eef2f4}footer{margin-top:35px;font-size:14px}@media(max-width:700px){main{margin:0;padding:20px}}"
    nav=''.join('<a href="'+RAW+n+'">'+label+'</a>' for n,label in [('CHAT_HANDOFF.md','给Chat入口'),('THEORY_DERIVATIONS.md','公式和证明'),('PAPER_DRAFT.md','论文初稿'),('BENCHMARKS.md','实验合同'),('COST_LEDGER.md','完整费用'),('REPRODUCE.md','复现')])
    page='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MAX-STRUCT-R7：低成本A+与有理材料搜索</title><style>'+css+'</style><main><nav>'+nav+'</nav>'+body+'<footer><a href="'+RELEASE+'">完整证据分卷</a> · <a href="../../">研究主页</a></footer></main></html>'
    missing=[u for u in re.findall(r'(?:src|href)="([^"#]+)"',page) if '://' not in u and not (ROOT/u).exists()]
    if missing:raise ValueError(missing)
    (ROOT/'index.html').write_text(page)
    (ROOT/'audit/HTML_EXPORT_RECEIPT.json').write_text(json.dumps(dict(status='EXPORTED',wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu,new_physical_calls=0,missing_local_links=missing),indent=2)+'\n')
    print(json.dumps(dict(status='EXPORTED',images=len(re.findall('<img ',page)),bytes=len(page.encode()))))
if __name__=='__main__':main()
