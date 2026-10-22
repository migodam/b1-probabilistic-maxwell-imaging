"""Separate offline gallery only: four predeclared DEV chi truth panels.

Does not rescore risk, alter frozen evaluation, run Target/physics/network or
train anything. Saved online arrays are gathered before these four truth reads.
"""
from pathlib import Path
import html
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import Book
from evaluate_p3l_joint_gaussian import PILOT,gather,figures,vector,write

OUT=ROOT/'research/p3l_local_posterior_imaging'
DEST=OUT/'B/joint_gaussian_figures'
JOB='plot_joint_gaussian_images'


def main():
    if (DEST/'GALLERY.json').exists():
        raise FileExistsError('gallery already recorded; no implicit plot replay')
    closure=json.loads((OUT/'B/jobs/joint_probit_gaussian/CLOSED.json').read_text())
    if closure.get('status')!='COMPLETE':
        raise RuntimeError('completed online joint Gaussian job required')
    risk_path=OUT/'B/JOINT_GAUSSIAN_IMAGE_RISK.json'
    if not risk_path.is_file():raise RuntimeError('dedicated risk evaluation must already exist')
    frozen_path=OUT/'configs/B_JOINT_GAUSSIAN_EVALUATION.json'
    frozen_text=frozen_path.read_text()
    frozen=json.loads(frozen_text)
    if frozen.get('IS_pilot')!=PILOT or PILOT!=['DEV_000','DEV_028','DEV_030','DEV_026']:
        raise RuntimeError('fixed pilot identities differ; no online/GT selection')
    declarations=[r for parent in PILOT for r in frozen['parents'] if r['parent']==parent]
    if [r['parent'] for r in declarations]!=PILOT:raise RuntimeError('fixed parent declaration missing')
    book=Book('B',JOB,max_wall=300,analysis_only=True)
    status='COMPLETE';error=None;started=time.perf_counter()
    try:
        online=json.loads((OUT/'B/JOINT_GAUSSIAN_REFERENCE.json').read_text())
        online_rows={r['parent']:r for r in online['rows']}
        estimates,rows=gather(OUT,declarations,online_rows)
        # Freeze which estimates/missing methods were present before GT reads.
        write(DEST/'GALLERY_INPUTS.json',dict(status='SAVED_ONLINE_ESTIMATES_LOADED_BEFORE_TRUTH',
            parents=PILOT,rows=rows,GT_read=False,GT_selection=False,analysis_only=True,
            figure_rule='shared real epsilon scale across all parents/methods; shared imaginary scale; separate uncalibrated modulus SD'))
        truth={}
        truth_reads=[]
        for parent in PILOT:
            path=OUT/f'inputs/offline_only/{parent}.npz'
            truth[parent]=vector(path,'chi')
            truth_reads.append(dict(parent=parent,path=path.relative_to(OUT).as_posix(),
                fields=['chi'],purpose='explicitly authorized offline gallery only',online_use=False))
        report=figures(OUT,estimates,truth)
        if frozen_path.read_text()!=frozen_text:raise RuntimeError('frozen evaluation changed during gallery; never overwritten here')
        limits=('原始四个 DEV 对象在看 GT 前固定；GT 仅用于本独立离线图示。'
            'GN Gaussian 图像与 IS 诊断均未取得完整 posterior 资格；uncalibrated modulus SD 不是校准后验不确定性或 credible interval。'
            '本脚本不重新评估风险，不调整估计器、quadrature、提案、parent 或冻结配置。')
        body='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P3-L · 固定 DEV 离线材料图库</title><style>body{font:17px/1.8 system-ui,sans-serif;background:#f6f8fa;color:#20323e;margin:0}main{max-width:1500px;margin:auto;padding:32px 20px}section{background:white;border:1px solid #dce4e8;border-radius:12px;padding:24px;margin:24px 0}img{width:100%;height:auto}a{color:#126c93}.notice{border-left:5px solid #b57334}code{background:#edf2f5}</style></head><body><main><nav><a href="../../index.html">P3-L 进展</a> · <a href="../JOINT_GAUSSIAN_IMAGE_RISK.json">已保存独立风险评价</a> · <a href="GALLERY.json">离线图库记录</a></nav><h1>固定四份 DEV 材料：联合 GN Gaussian 与强基线</h1><section class="notice"><p>'''+html.escape(limits)+'''</p><p>行顺序：DEV_000 / DEV_028 / DEV_030 / DEV_026。列顺序：离线真值、历史 native-q anchor、probit 坐标 mode、该 mode 固定属性 Gaussian mean、joint GH17 mean、未合格 full-target IS mean。缺失结果明确显示 MISSING / FAILED，不替换样本。</p><p>Re ε = 1 + Re χ；Im ε = Im χ。两幅材料图分别在全部 parent / method 间使用统一色标。色标只服务离线显示，不反馈到任何在线计算。</p></section>'''
        for key,title in [('real','Re epsilon：统一色标'),('imag','Im epsilon：统一色标'),('uncertainty','GN Gaussian modulus SD：未校准')]:
            item=report[key]
            if 'path' not in item:raise RuntimeError(f'gallery figure missing: {key}')
            filename=Path(item['path']).name
            body+='<section><h2>'+html.escape(title)+'</h2><img src="'+html.escape(filename,quote=True)+'" alt="'+html.escape(title,quote=True)+'"><p>显示范围 ['+str(item['shared_min'])+', '+str(item['shared_max'])+']。</p>'
            if key=='uncertainty':body+='<p>sqrt(E|χ−Eχ|²)，是给定 GN Gaussian reference 的复模标准差；不是 Re/Im 分量标准差，也不代表经过验证的 posterior spread。</p>'
            body+='</section>'
        body+='<section><h2>证据边界</h2><p>独立 plot job 只支付纯离线分析 wall；新增物理 RHS、state、LU、GPU 与训练均为 0。材料真值读取限于声明四个 development parent 的 chi，未读取 locked 数据。</p><p><a href="GALLERY_INPUTS.json">GT 读取前的估计器清单</a> · <a href="../jobs/'+JOB+'/CLOSED.json">分析作业闭合收据</a> · <a href="../../configs/B_JOINT_GAUSSIAN_EVALUATION.json">原冻结评价配置</a></p></section></main></body></html>'
        (DEST/'index.html').write_text(body)
        receipt=dict(status='COMPLETE_SEPARATE_OFFLINE_GALLERY_ONLY',job=JOB,parents=PILOT,
            figures=report,truth_reads=truth_reads,GT_selection=False,online_GT_use=False,
            risks_recomputed=False,frozen_evaluation_modified=False,online_estimates_loaded_before_truth=True,
            shared_color_scales=True,uncertainty_qualification='UNCALIBRATED_GAUSSIAN_REFERENCE_NOT_POSTERIOR',
            new_physical_RHS=0,new_states=0,new_LU=0,GPU_seconds=0,training_runs=0,
            analysis_only=True,gallery_wall_seconds=time.perf_counter()-started,limitations=limits)
        write(DEST/'GALLERY.json',receipt)
        print(json.dumps(dict(status=receipt['status'],gallery=str(DEST/'index.html'),
            figures=[str(OUT/v['path']) for v in report.values()],new_physical_RHS=0)))
    except BaseException as exc:
        status='ERROR';error=repr(exc)
        write(DEST/'GALLERY_FAILURE.json',dict(error=error,new_physical_RHS=0,analysis_only=True))
        raise
    finally:
        book.close(status,error)


if __name__=='__main__':main()
