"""Summarize already sealed/scored historical rows; opens no truth files."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main():
    out=ROOT/'RESULTS/history_warm'
    scored=json.loads((out/'SCORES.json').read_text())
    assert scored['status']=='HISTORICAL_ONLY_SCORED_AFTER_ALL21_SEALED'
    rows=scored['rows'];ids=list(dict.fromkeys(r['case_id'] for r in rows))
    methods=['R8_nn_cie_dual','R8_nn_a_plus','R8_nn_a_plus_sourcewise','R8_nn_physical_lbfgs','R8_uniform_a_plus_sourcewise']
    summary=[]
    for method in methods:
        r=[r for r in rows if r['method']==method]
        summary.append(dict(method=method,objects=len(r),historical_material_successes=sum(x['material_image_success'] for x in r),CPU_seconds=sum(x['cpu_seconds'] for x in r),LU=sum(x['LU'] for x in r),complex_RHS=sum(x['complex_RHS'] for x in r)))
    (out/'PAIRED_SUMMARY.json').write_text(json.dumps(dict(scope='FOUR_EXPOSED_HISTORICAL_OBJECTS_ONLY',source='SCORES.json',new_truth_or_physics_calls=0,method_counts=summary,rows=[dict(case_id=r['case_id'],method=r['method'],amplitudeNRMSE=r['material_nrmse'],IoU=r['otsu_iou'],original_Phi=r['original_Phi'],stop=r['corrector_stop'],iterations=r['corrector_iterations'],CPU=r['cpu_seconds']) for r in rows if r['method'].startswith('R8_')]),indent=2))
    lines=['','## 有限历史回归的结果解释','',
      '这四个旧对象上，冻结原10通道、材料初始化且APCS初始为零的NN没有替代A更新在两个困难弯条上的收益：NN+Dual为0/4成功；NN+global A+、NN+sourcewise A+各2/4；NN+direct physics为1/4；uniform sourcewise A+也是2/4。这里成功仅指原历史材料指标，不是新确认、统计泛化或唯一性结论。','',
      '- 7b8d：旧uniform original A和R7 final global A+成功，新增uniform sourcewise A+也成功。NN+global/sourcewise A+保留成功（NRMSE .04559/.04551，IoU均.94382），NN+Dual仍停在错误结构（.39310/.40777），NN+direct也失败（.40726/.23762）。NN没有把Dual移到A达到的材料解。','- b20d：NN+global/sourcewise A+保留成功（.04096/.04179，IoU .95181/.94048），NN+Dual失败（.56656/.25000）。NN+direct独立成功（.04102/.87912），这是本回归中唯一不使用A的NN成功；不能将它推广到另外三个对象。','- 42cb：同版本同预算新增uniform global A+为.75888/.11111，NN+global A+为.70839/.13158，两者均失败。NN+Dual、NN+sourcewise A+和NN+direct亦全部失败；NN direct的较低NRMSE .43324伴随IoU仅.16964和residual .40294，不能称为恢复。','- 40ac开缝壳：所有旧基线及本轮四NN臂均失败。NN+global A+ .62595/.27660与R7 uniform global A+ .62597/.27660几乎相同；NN+sourcewise .62416/.27128没有修复拓扑。','',
      '四个NN初始材料NRMSE均比uniform更低，初始原Phi也均较低；最终42cb与40ac仍失败，说明本回归中的初始接近并不足以保证原非线性校正恢复。所有初值均完整执行，未因初始Phi拒绝。','',
      '源方向选择的效应也有条件：sourcewise A+没有在这四对象上增加成功数；其CPU/终止路径随对象和初始化不同。逐方法的原Phi、corrector实际步数、终止原因和CPU已在完整表中列出。PGTOL/relative-F是数值停止，并非材料成功证书；42cb部分路径达到3000步上限，不称为充分收敛。','',
      '比较采用相同预设最大预算，实际CPU并不相等。R7旧uniform global A+在两个bar上分别45.36/46.53 CPU秒；NN global A+为57.23/50.89秒，此处没有测得NN加速。新增uniform sourcewise A+两个bar为53.17/65.28秒，NN sourcewise为52.45/57.31秒；这些有限路径差异不足以作通用加速结论。direct最多6000步，其停止与辅助臂3000+3000预算不同。','',
      '全批21项solver CPU1274.461107秒、任务wall之和1583.460057秒，串行wrapper CPU1275.020020秒、wall1584.147680秒；前处理另0.077402 CPU秒、0.142219 wall秒。计入LU40204、复数RHS964056；共享features含4次receiver SVD，公开Green动作/接收adjoint各48源列，无额外LU。原网络训练387.690002 CPU秒及训练数据生成费用属于复用外部成本，不藏入零成本。','',
      '上述结果只检验固定的material-only warm start。未来保留NN物理APCS偏置的初始化属于另一协议；本次没有补跑、重选或根据历史GT改动当前21项。','']
    p=ROOT/'HISTORICAL_WARM_STARTS.md';text=p.read_text();marker='\n## 有限历史回归的结果解释'
    if marker in text:text=text.split(marker)[0]
    p.write_text(text+'\n'.join(lines))
    (out/'FIGURE_REVIEW.json').write_text(json.dumps(dict(status='VISUALLY_REVIEWED',figure='FIGURES/history_warm/all_four_objects.png',objects=4,columns=10,all_four_objects_and_prespecified_method_columns_present=True,full_common_row_material_range=True,failed_spikes_to20_visible=True,new_physical_calls=0),indent=2))
    print(json.dumps(dict(status='INTERPRETED_SEALED_HISTORICAL_ROWS',report=str(p),methods=summary)))

if __name__=='__main__':main()
