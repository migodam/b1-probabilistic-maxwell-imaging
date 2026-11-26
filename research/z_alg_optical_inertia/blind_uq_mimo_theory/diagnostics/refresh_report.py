"""Fill the two report tables using only scored, independently replayed results."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / 'results'

def csv_rows(name):
    with (out / name).open() as stream:
        return list(csv.DictReader(stream))

def between_markers(doc, label, body):
    start = '<!-- LIVE-' + label + ' -->'
    end = '<!-- END-LIVE-' + label + ' -->'
    first, section = doc.split(start)
    _, last = section.split(end)
    return first + start + '\n' + body + '\n' + end + last

regions = csv_rows('region_intervals.csv')
cases = csv_rows('objects.csv')
summary = json.loads((out / 'region_summary.json').read_text())
s = summary['all']
top = (
    f"在主 8 发、8 收读数下，**{s['exact_hull']}/{s['regions']} 个区域最窄的保证区间恰为 [0, 4.5]**。"
    f"32² 是 {summary['n32']['exact_hull']}/{summary['n32']['regions']}，原生 96² 是 {summary['n96']['exact_hull']}/{summary['n96']['regions']}。"
    f"其中 {s['continuous']} 个区域还证明了从 0 到 4.5 的每一个峰值都实际可取。"
    f"解封真值后发现，{s['removable_with_real_cells']} 个可以归零的亮区里本来都有材料。"
    "因此，一张归零后仍能解释测量的图可以作为替代解释，不能作为实际样品内该处为空气的判据。"
)
table = [
    '| 对象 | 盲选区域数 | 两端均可达 | 整段均可达 | CIE 更新次数 | 原全波残差 / 半径 | 全场峰值严格大于 |',
    '|---|---:|---:|---:|---:|---:|---:|',
]
for c in cases:
    rows = [r for r in regions if r['case'] == c['case']]
    endpoints = sum(r['width_exact'] == 'True' for r in rows)
    paths = sum(r['full_continuous'] == 'True' for r in rows)
    table.append(f"| {c['case']} | {len(rows)} | {endpoints} | {paths} | {int(c['baseline_updates']):,} | {float(c['stationary_residual_over_r']):.2f} | {float(c['global_peak_lower']):g} |")
target = ROOT / 'BENCHMARK_REPORT_ZH.md'
article = between_markers(target.read_text(), 'COUNTS', top)
target.write_text(between_markers(article, 'OBJECTS', '\n'.join(table)))
