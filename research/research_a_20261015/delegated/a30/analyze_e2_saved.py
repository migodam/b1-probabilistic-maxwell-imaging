"""Saved-array-only E2 analysis. Does not import any physics implementation."""
from time import perf_counter, process_time
_start, _cpu = perf_counter(), process_time()
from pathlib import Path
import json
from collections import Counter
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "research/research_a_20261015/results/E2"
STATUS = "FAILED"
PREVIOUS_RECEIPT = json.loads((OUT/"E2_POSTPROCESS_JOB.json").read_text()) if (OUT/"E2_POSTPROCESS_JOB.json").exists() else None


def write(name, data):
    (OUT/name).write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def mean(values):
    return float(np.mean(values)) if values else None


def mdistance(a,b):
    return float(np.sqrt(np.mean(np.abs(a-b)**2))/1.2)


try:
    history = json.loads((OUT/"SEARCH_HISTORY.json").read_text())
    metrics = json.loads((OUT/"METRICS.json").read_text())
    protocol = json.loads((OUT/"PROTOCOL.json").read_text())
    result = json.loads((OUT/"RESULTS.json").read_text())
    ledger = [json.loads(s) for s in (OUT/"COST_LEDGER.jsonl").read_text().splitlines()]
    audit = json.loads((OUT/"ACTUAL_BACKEND_AUDIT.json").read_text())
    old_failure = json.loads((OUT/"RESULTS_FAILED_AUDIT_ATTEMPT1.json").read_text())
    revision = json.loads((OUT/"AUDIT_REVISION.json").read_text())
    cases, arms = protocol["cases"], protocol["arms"]
    mm = {(r["case_id"],r["arm"]):r for r in metrics}
    threshold = protocol["compatibility_radius"]
    case_results=[]
    trajectories=[]
    added=[]
    with np.load(OUT/"MATERIAL_AND_PREDICTION_ARRAYS.npz", allow_pickle=False) as arrays:
        for case in cases:
            base=arrays[case+"__base"]
            bg=arrays[case+"__background"]
            truth=arrays[case+"__OFFLINE_TRUTH"]
            base_mask=np.abs(base-bg)>.18
            base_row=next(r for r in history if r["case_id"]==case and r["arm"]=="base")
            perarm={}
            for arm in arms:
                metric=mm[case,arm]
                rows=[r for r in history if r["case_id"]==case and r["arm"]==arm]
                archive=[]
                if base_row["compatible"]:
                    archive.append({"name":"base","chi":base,"residual":base_row["residual"]})
                accepted_trials=0
                compatible_trial_scores=[]
                for row in rows:
                    if not row["compatible"]:
                        continue
                    accepted_trials+=1
                    image=arrays[row["array_prefix"]+"__chi"]
                    compatible_trial_scores.append(dict(name=row["array_prefix"],
                        full_complex_MAE=float(np.mean(np.abs(image-truth))),
                        distance_to_base=mdistance(image,base),
                        structural_mask_changed_cells=int(np.count_nonzero((np.abs(image-bg)>.18)!=base_mask))))
                    if all(mdistance(image,r["chi"])>=.03 for r in archive):
                        archive.append({"name":row["array_prefix"],"chi":image,"residual":row["residual"]})
                if len(archive)!=metric["valid_distinct_material_count"]:
                    raise AssertionError("Archive reconstruction disagrees with saved metric")
                details=[]
                for item in archive:
                    image=item["chi"]
                    mask=np.abs(image-bg)>.18
                    details.append(dict(name=item["name"],residual=item["residual"],
                        material_distance_to_base=mdistance(image,base),
                        changed_structure_mask_cells=int(np.count_nonzero(mask!=base_mask)),
                        full_complex_MAE=float(np.mean(np.abs(image-truth))),
                        full_Re_MAE=float(np.mean(np.abs((image-truth).real))),
                        full_Im_MAE=float(np.mean(np.abs((image-truth).imag))),
                        relative_complex_MAE_improvement_over_base=float(1-np.mean(np.abs(image-truth))/np.mean(np.abs(base-truth))),
                        structural_minimum6cells_met=bool(np.count_nonzero(mask!=base_mask)>=6)))
                    if item["name"]!="base":
                        added.append(dict(case_id=case,arm=arm,**details[-1]))
                trajstats=[]
                for t in sorted({r["trajectory"] for r in rows}):
                    tr=sorted([r for r in rows if r["trajectory"]==t],key=lambda r:r["step"])
                    first,last=tr[0],tr[-1]
                    first_image=arrays[first["array_prefix"]+"__chi"]
                    last_image=arrays[last["array_prefix"]+"__chi"]
                    stat=dict(case_id=case,arm=arm,trajectory=t,round=first.get("round"),
                        real_parameters=first.get("real_parameters"),
                        first_residual=first["residual"],last_residual=last["residual"],
                        best_residual=min(r["residual"] for r in tr),
                        residual_reduction_fraction=float(1-last["residual"]/first["residual"]),
                        initial_compatible=first["compatible"],any_compatible=any(r["compatible"] for r in tr),
                        invalid_to_valid=not first["compatible"] and any(r["compatible"] for r in tr),
                        material_correction_distance=mdistance(first_image,last_image),
                        protected_cells=first.get("protected_cells"),
                        protected_defect_max=max(r.get("protected_defect",0.) for r in tr),
                        box_projection_defect=first.get("edit",{}).get("box_projection_defect"))
                    trajectories.append(stat)
                    trajstats.append(stat)
                perarm[arm]=dict(fresh_states=len(rows),compatible_full_trials=accepted_trials,
                    valid_distinct_material_count=len(archive),
                    structural_families=metric["structurally_distinct_masks_at_least6changedcells"],
                    best4_complex_MAE=metric["Best_of_4_complex_MAE"],
                    selected_complex_MAE=(metric["selected_offline_metrics"] or {}).get("full_complex_MAE"),
                    coverage4=metric["Coverage_at_4"],archive=details,
                    all_compatible_trial_offline_scores=compatible_trial_scores,
                    pre_dedup_offline_oracle_complex_MAE=min([r["full_complex_MAE"] for r in compatible_trial_scores]+
                        ([float(np.mean(np.abs(base-truth)))] if base_row["compatible"] else []),default=None),
                    invalid_to_valid_trajectories=sum(t["invalid_to_valid"] for t in trajstats),
                    best_trial_residual=min(r["residual"] for r in rows))
            case_results.append(dict(case_id=case,parent_id=mm[case,arms[0]]["parent_id"],
                family=mm[case,arms[0]]["family"],base_residual=base_row["residual"],base_compatible=base_row["compatible"],
                base_full_complex_MAE=float(np.mean(np.abs(base-truth))),
                arms=perarm))
    parent_results=[]
    for parent in dict.fromkeys(x["parent_id"] for x in case_results):
        pair=[r for r in case_results if r["parent_id"]==parent]
        parent_results.append(dict(parent_id=parent,family=pair[0]["family"],conditions=len(pair),
            base_compatible_conditions=sum(r["base_compatible"] for r in pair),
            base_full_complex_MAE=mean([r["base_full_complex_MAE"] for r in pair]),
            arms={arm:dict(mean_valid_distinct_material_count=mean([r["arms"][arm]["valid_distinct_material_count"] for r in pair]),
                mean_structural_families=mean([r["arms"][arm]["structural_families"] for r in pair]),
                compatible_conditions=sum(r["arms"][arm]["valid_distinct_material_count"]>0 for r in pair),
                coverage4_conditions=sum(r["arms"][arm]["coverage4"] for r in pair),
                paired_conditions_Best4_complete=all(r["arms"][arm]["best4_complex_MAE"] is not None for r in pair),
                paired_conditions_Best4_mean=mean([r["arms"][arm]["best4_complex_MAE"] for r in pair])
                    if all(r["arms"][arm]["best4_complex_MAE"] is not None for r in pair) else None)
                for arm in arms}))
    summaries={}
    for arm in arms:
        tr=[r for r in trajectories if r["arm"]==arm]
        am=[r["arms"][arm] for r in case_results]
        summaries[arm]=dict(full_trial_states=sum(r["fresh_states"] for r in am),
            compatible_full_trials=sum(r["compatible_full_trials"] for r in am),
            valid_distinct_material_count_by_case=[r["valid_distinct_material_count"] for r in am],
            structural_families_by_case=[r["structural_families"] for r in am],
            total_archive_members_including_shared_base=sum(r["valid_distinct_material_count"] for r in am),
            new_archive_members_beyond_shared_base=sum(r["valid_distinct_material_count"]-int(c["base_compatible"]) for r,c in zip(am,case_results)),
            any_valid_cases=sum(r["valid_distinct_material_count"]>0 for r in am),
            no_valid_cases=sum(r["valid_distinct_material_count"]==0 for r in am),
            best4_all8conditions_mean=None,
            best4_eligible3conditions_descriptive_mean=mean([r["best4_complex_MAE"] for r in am if r["best4_complex_MAE"] is not None]),
            coverage4_cases=sum(r["coverage4"] for r in am),
            invalid_to_valid_trajectories=sum(r["invalid_to_valid"] for r in tr),
            total_trajectories=len(tr),
            median_initial_residual=float(np.median([r["first_residual"] for r in tr])),
            median_final_residual=float(np.median([r["last_residual"] for r in tr])),
            median_paired_residual_reduction_fraction=float(np.median([r["residual_reduction_fraction"] for r in tr])),
            median_material_correction_distance=float(np.median([r["material_correction_distance"] for r in tr])),
            max_protected_defect=max(r["protected_defect_max"] for r in tr),
            median_box_projection_defect= float(np.median([r["box_projection_defect"] for r in tr if r["box_projection_defect"] is not None]))
                if any(r["box_projection_defect"] is not None for r in tr) else None)
    qps=[r for r in ledger if r["event"]=="MATERIAL_QP"]
    executions=[r for r in ledger if r["event"]=="EXECUTION_END"]
    costs=dict(result["resources"],remaining_RHS=9216-result["resources"]["new_complex_LU_RHS"],
        remaining_states=384-result["resources"]["new_full_state_attempts"],new_physics_this_postanalysis=0,
        new_RHS_this_postanalysis=0,failed_audit_attempts=old_failure["resources"]["new_full_state_attempts"],
        failed_audit_RHS=old_failure["resources"]["new_complex_LU_RHS"],
        execution_wall_from_append_only_ledger=sum(r["execution_wall_seconds"] for r in executions))
    output=dict(status="SAVED_ARRAY_ANALYSIS_COMPLETE",scientific_scope="8exposed conditions/4parents;2DscalarTM",
        threshold=threshold,costs=costs,backend_audit=audit,
        QP_diagnostics=dict(total=len(qps),converged=sum(r["converged"] for r in qps),iteration_caps=sum(not r["converged"] for r in qps)),
        arms=summaries,cases=case_results,parents=parent_results,new_archive_members=added,
        trajectories=trajectories,
        conclusions=dict(new_structural_family_beyond_multistart=0,
            compensation_created_new_valid_finite_edit_trajectories=summaries["finite_history"]["invalid_to_valid_trajectories"],
            best4_improvement_on_all_eligible_cases=0.,all8_best4_mean="UNDEFINED5/8NO_VALID_SOLUTION",
            physical_global_absence="NOT_ESTABLISHED",learned_A3="NOT_RUN",E3="NOT_AUTHORIZED",
            next_action="Stop this compensation superiority claim; fix coarse reconstruction/support/longer correction only under new declaredbudget"),
        protocol_deviation=revision,
        no_probability_claim=True,offline_saved_truth_read=True,online_truth_used=False,locked_payload_read=False)
    write("RESULT_SUMMARY.json",output)
    table=[]
    for r in case_results:
        aa=r["arms"]
        table.append(f"|{r['case_id']}|{r['base_residual']:.4f}|{aa['finite_history']['valid_distinct_material_count']} / {aa['multistart']['valid_distinct_material_count']} / {aa['weak_null']['valid_distinct_material_count']} / {aa['same_edit_no_compensation']['valid_distinct_material_count']}|{aa['finite_history']['structural_families']}|"+
                     (f"{aa['finite_history']['best4_complex_MAE']:.8f}" if aa['finite_history']['best4_complex_MAE'] is not None else "NOT_DEFINED")+"|")
    extra=added[0] if added else None
    finite=summaries["finite_history"]
    descriptive=finite["best4_eligible3conditions_descriptive_mean"]
    text=f'''# E2：已保存 full-Maxwell 数组的失败分析

本文件只分析已付费数组，不调用 solver、不新增 RHS，不更改原 E2、A30 或失败结果。统计单位为4个 parent，各 weak/ring 两条件。全部是历史已暴露 development；不是独立验证。

## 主要结果

**没有找到比普通 multi-start 更多的结构解族；有限编辑补偿没有把任何原先不相容的编辑修正成新相容轨迹。Best-of-4 在所有具有合法候选的三个条件上与普通 multi-start、weak/null、无补偿均相同；其余5/8条件没有合法解，不能删掉后报告总体成功。**

有限编辑 archive 计数为`{finite['valid_distinct_material_count_by_case']}`；multi-start 和 weak/null 都为`{summaries['multistart']['valid_distinct_material_count_by_case']}`。有限编辑多出的一个是**材料场 RMS 去重**意义的额外成员；同样成员在无补偿的原始编辑中已经存在。它不满足预注册的至少6个结构 mask 改变单元，不能称新结构族，更不能归功于 nonlinear compensation。

|条件|base白化残差|合法去重材料数 finite / multistart / weak / no-comp|finite结构族数|三方法及no-comp相同的Best4复材料MAE|
|---|---:|---|---:|---:|
{chr(10).join(table)}

相容阈值固定为sqrt(χ²_512,.95)={threshold:.9f}；原 sqrt(2)/σ、8源32接收器、固定gain1、材料box均未放宽。base只有3/8满足相容性，所有ring条件都失败；finite-history只在1/8具有第二个去重材料，结构族仍均不超过1。四个parent没有一个在两个背景条件都具有合法候选，所以全parent双条件 Best4 均为NOT_DEFINED。不得将仅weak成功的三个parent写成完整paired endpoint。

## 额外材料到底是什么

唯一额外有限编辑 archive 成员见`RESULT_SUMMARY.json/new_archive_members`。它来自 p250069_weak 的第三条轨迹，初始 full residual 已经低于相容阈值，无需补偿。该材料与base的 field-RMS距离达到.03去重标准，但阈值结构 mask 改变小于6格。其复材料MAE比base更差；TopK离线择优仍选择base。于是“多一个材料成员”不等于“多发现一个重要结构”，也不等于“成像改善”。no-comp保留同一成员，是直接消融反证。

三个有合法候选条件的 Best4 描述性均值为{descriptive:.10f}，各方法相同，改善0。这个数仅描述三个已成功条件；5个无合法解的MAE没有用0、无穷或base偷偷填补。Coverage@4按原离线field-distance≤.1定义，见完整表及JSON；该coverage门槛与材料去重.03、结构6格是不同定义。

必须保留一个去重策略诊断：如果离线查看**全部合法trial、尚未按.03去重**，multi-start在p250069_weak的一个trial使MAE从.04351438到.04308716，改善.9818%；finite-history在p250077_weak的一个原始编辑从.04845289到.04803630，改善.8598%。它们被在线first-member材料去重并入base邻域，因此没有改变已冻结archive的Best4。这些微小的within-family改善不构成新结构或finite-compensation优势；前者反而提醒不能把“archive无改善”误说成“所有物理trial的材料完全没有改善”。未经预注册的truth择优只是事后诊断，不能回写在线筛选规则。

## 修正器确实改变了材料，但未完成目标

共{len(qps)}个小材料QP，{sum(r['converged'] for r in qps)}个报告收敛，{sum(not r['converged'] for r in qps)}个达到迭代上限。full-forward修正产生真实非零材料变化，并非接口空转。有限编辑32条轨迹的初始/最终 residual 中位数分别为{finite['median_initial_residual']:.4f}/{finite['median_final_residual']:.4f}，配对残差下降中位比例{finite['median_paired_residual_reduction_fraction']:.2%}，但 invalid→valid 轨迹数为0/32。protected edit cells 最大缺陷{finite['max_protected_defect']:.3g}，优化器没有偷偷把洞填回去。

普通 multi-start 和 weak/null 的全局起点更远，经过两个修正步骤后大幅降低 residual；部分进入相容阈值，但都落在base的.03材料邻域内。它们没有发现额外结构。这支持“当前局部修正/表示预算不足以产生可验证新分支”的解释；不能从中推导这些方法的全局最佳能力。weak/null使用的是N3处材料tangent，且有限步还经过材料box投影，投影与有限非线性都会破坏线性弱方向性质；它不是任意receiver-dark current，更不是已经验证的有限材料null族。

Round1为16实参数低频DCT，Round2为32实参数multiscale Haar；只运行每轨迹两次 nonlinear correction。材料是576格的完整χ，源共享；chart不是一般高维材料空间的完备表示。A30旧sparse64的35.69%仍是truth-fit容量，不能套用到本轮16/32在线优化。五个无解条件的base本身不相容，说明本轮连共同粗重建前置条件也未充分满足。缺少可行base、有限参数表示、有限迭代和物理不可补偿这几种原因，在此预算内没有被完全分离。

## parent汇总与成本

`RESULT_SUMMARY.json/parents`先平均每parent的weak/ring。finite合法材料数的parent均值依次0/1/.5/.5；对照0/.5/.5/.5；finite唯一+.5parent差对应上文无结构收益的同一个原始编辑。结构族parent均值各方法一致0/.5/.5/.5。全8条件算法新增结构族数0；与multi-start相比新增0；best4收益0或未定义。没有用32条轨迹当32个独立parent，也没有为4个暴露parent编造确认性置信区间。

累计{costs['new_full_state_attempts']} /384新full-state attempts、{costs['new_complex_LU_RHS']} /9216complex LU RHS，包含首次失败audit的1state/91RHS；受控执行wall {costs['controlled_execution_wall_seconds']:.5f}s。剩56state但只有2RHS，无法支持下一个8源full state。无需继续付费探针。GPU0，旧P2额度0，旧A30 ledger不改。

成功audit的实际finite secant、Schur、dual、材料mass-adjoint和c=χe误差均<1e−8，见ACTUAL_BACKEND_AUDIT.json。首次dual小二阶量相对归一化误差1.91668e−8失败保留；parent改为可解析1%扰动的dual试验，容差不变，不把失败抹掉。为保留91RHS重复费用，省略两个原拟定no-comp端点（最初两个条件各少一次）；完整运行no-comp30states而不是32。原PROTOCOL仍保留原计划，实际偏差以AUDIT_REVISION和ledger为准。五方法/全预算speed优越性没有建立。

本分析脚本另记`E2_POSTPROCESS_JOB.json`，0physics/0RHS；其耗时不混入已冻结的18.13s求解stage。

## 科学判断边界

本轮否定的是**已实施低维受限corrector和编辑组合在该8条件下的新增结构/Best4优势**。它没有证明完整合法材料空间不存在其他相容结构，没有证明A30所有局部法一定失败，也没有证明多源TM逆问题全局唯一。没有概率认证；candidate-only字段应为null。A3 learned policy未训练，有限history parent重用未证明收益。没有比较TSOM/ROM、反馈特化或learned correction，不能为它们归因。

E3未授权；这轮也没有提供由二维成功机制进入三维的科学理由。优先停止当前补偿优越性主张，保留失败作为反证。若未来另获预算，应先解决强粗重建与真正可行base，再用同预算独立起点和扩展表示分离优化瓶颈；不能靠扩noise/freegain或加大网络制造多解。最终机制裁决由parent Codex负责。
'''
    (OUT/"E2_FAILURE_ANALYSIS.md").write_text(text)
    STATUS="COMPLETE"
finally:
    current_run=dict(status=STATUS,wall_seconds=perf_counter()-_start,process_cpu_seconds=process_time()-_cpu)
    previous_runs=(PREVIOUS_RECEIPT.get("runs",[dict(status=PREVIOUS_RECEIPT["status"],wall_seconds=PREVIOUS_RECEIPT["wall_seconds"],
        process_cpu_seconds=PREVIOUS_RECEIPT["process_cpu_seconds"])]) if PREVIOUS_RECEIPT else [])
    runs=previous_runs+[current_run]
    write("E2_POSTPROCESS_JOB.json",dict(status=STATUS,role="saved_array_postanalysis",new_physics_attempts=0,
        new_RHS=0,GPU_seconds=0,old_ledger_modified=False,wall_seconds=perf_counter()-_start,
        process_cpu_seconds=process_time()-_cpu,runs=runs,
        cumulative_postprocess_wall_seconds=sum(r["wall_seconds"] for r in runs),
        cumulative_postprocess_cpu_seconds=sum(r["process_cpu_seconds"] for r in runs),
        scope="This Python saved-array statistics/document generation only;interactive review excluded",
        outputs=["E2_FAILURE_ANALYSIS.md","RESULT_SUMMARY.json"],
        inputs=["SEARCH_HISTORY.json","METRICS.json","MATERIAL_AND_PREDICTION_ARRAYS.npz","COST_LEDGER.jsonl","ACTUAL_BACKEND_AUDIT.json"]))

print(json.dumps(dict(status=STATUS,costs=costs,arms={k:{x:v for x,v in a.items() if x in ["valid_distinct_material_count_by_case","structural_families_by_case","invalid_to_valid_trajectories","coverage4_cases"]} for k,a in summaries.items()},extras=added)))
