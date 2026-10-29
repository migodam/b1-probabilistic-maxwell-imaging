"""E3 zero-new-PDE old-bank conditional consensus, then sealed offline scoring."""
from pathlib import Path
import argparse
import csv
import io
import json
import os
import resource
import struct
import sys
import time
import zipfile
import zlib

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT.parent / "research_w_gq"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(OLD / "src"))
sys.path.insert(0, str(OLD / "tools"))
from p3_bank import read_online_prior, load_keys, safe_path, PRIOR_ID
from r2.consensus import (LABEL, METHODS, SIGMA, TV_GATE, make_schedule, online_core,
                          algebra_checks, conditional_reference, static_distances, offline_distances)


def plain(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, dict):
        return {str(k): plain(v) for k,v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    return value


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(plain(value), stream, indent=2, allow_nan=False)
        stream.write("\n")


def write_csv(path, rows):
    with Path(path).open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def event(book, label, **meta):
    if book is not None:
        book.event("E3", label, **meta)


class BaseRowOracle:
    """Read only the leading TX0/1 bytes of data.npy inside each paid row.

    The legacy row member stores TX0..3.  Loading its entire ndarray would
    expose outer TX2/3, so this provider parses the NPY header and streams only
    the leading two TX rows.  No outer array values are loaded or interpreted.
    """
    def __init__(self, root):
        self.directory = safe_path(root) / "paid_rows"
        self.calls = []
        if not self.directory.is_dir():
            raise FileNotFoundError("individual paid_rows are required; no full-bank fallback")

    def __call__(self, atom):
        if not isinstance(atom, (int,np.integer)) or not 0 <= atom < 2950:
            raise ValueError("legal original finite atom index required")
        path = self.directory / f"{int(atom):04d}.npz"
        with zipfile.ZipFile(path) as archive:
            with archive.open("row_index.npy") as stream:
                row = np.lib.format.read_array(stream, allow_pickle=False)
            if int(row) != int(atom):
                raise ValueError("paid-row atom identity mismatch")
            # ZipExtFile buffers decompressed bytes beyond read(n).  Bounded
            # raw-deflate decoding prevents that buffer exposing TX2/3 values.
            first = self.member_prefix(archive, "data.npy", 12)
            if first[:6] != b"\x93NUMPY":
                raise ValueError("NPY data member required")
            version=tuple(first[6:8])
            prefix=10 if version==(1,0) else 12 if version==(2,0) else 0
            if not prefix:
                raise ValueError("unsupported paid-row NPY schema")
            header_size=struct.unpack("<H" if prefix==10 else "<I",first[8:prefix])[0]
            header_end=prefix+header_size
            stream=io.BytesIO(self.member_prefix(archive,"data.npy",header_end))
            np.lib.format.read_magic(stream)
            reader=np.lib.format.read_array_header_1_0 if prefix==10 else np.lib.format.read_array_header_2_0
            shape,fortran,dtype=reader(stream)
            if shape!=(4,16) or fortran or dtype.kind!="c" or dtype.hasobject:
                raise ValueError("frozen C-order 4TX by16RX complex row required")
            size=2*16*dtype.itemsize
            payload=self.member_prefix(archive,"data.npy",header_end+size)[header_end:]
            if len(payload)!=size:
                raise ValueError("truncated base-port paid row")
            value=np.frombuffer(payload,dtype=dtype).reshape(2,16).copy()
        self.calls.append(int(atom))
        return value

    @staticmethod
    def member_prefix(archive,name,size):
        info=archive.getinfo(name); stream=archive.fp
        stream.seek(info.header_offset)
        header=stream.read(30)
        if len(header)!=30 or header[:4]!=b"PK\x03\x04":
            raise ValueError("ZIP local header required")
        filename_size,extra_size=struct.unpack("<HH",header[26:30])
        stream.seek(filename_size+extra_size,1)
        if info.compress_type==zipfile.ZIP_STORED:
            return stream.read(size)
        if info.compress_type!=zipfile.ZIP_DEFLATED:
            raise ValueError("unsupported paid-row compression")
        encoded=stream.read(info.compress_size)
        decoder=zlib.decompressobj(-zlib.MAX_WBITS)
        return decoder.decompress(encoded,max_length=size)


def save_online(directory, result, prior):
    directory.mkdir()
    arrays = {k:v for k,v in result.items() if isinstance(v,np.ndarray)}
    arrays.update(pose_nodes=prior.pose_nodes, wall_group_ids=prior.group_ids)
    np.savez_compressed(directory / "QUERIES_AND_CONDITIONALS.npz", **arrays)
    metadata = {k:v for k,v in result.items() if not isinstance(v,np.ndarray)
                and k not in ("methods","intervals","query_trace")}
    metadata["method_metadata"] = {}
    for method in METHODS:
        values = result["methods"][method]
        np.savez_compressed(directory / f"{method}.npz", **{k:v for k,v in values.items() if isinstance(v,np.ndarray)})
        metadata["method_metadata"][method] = {k:v for k,v in values.items() if not isinstance(v,np.ndarray)}
    interval = result["intervals"]
    np.savez_compressed(directory / "DEFAULT_UNKNOWN_TAIL_INTERVALS.npz", **{k:v for k,v in interval.items() if isinstance(v,np.ndarray)})
    metadata["default_unknown_tail_events"] = interval["events"]
    write_json(directory / "QUERY_TRACE.json", result["query_trace"])
    write_json(directory / "ONLINE.json", metadata)
    write_json(directory / "OUTPUTS_SEALED.json", dict(method_ids=METHODS,
        all_method_outputs_written=True, reference_or_truth_seen=False, outer_TX_seen=False,
        shared_unique_atom_queries=len(result["query_indices"]), label=LABEL))


def verify_seal(online):
    seal = json.loads((online / "ALL_OUTPUTS_SEALED.json").read_text())
    if (seal.get("case_count") != 24 or seal.get("reference_or_truth_seen") is not False
            or seal.get("all_method_outputs_written") is not True or tuple(seal.get("method_ids",())) != METHODS):
        raise ValueError("all 24 cases and three online methods must seal before scoring")
    for parent in range(12):
        for repeat in range(2):
            directory = online / f"parent_{parent:02d}_repeat_{repeat:02d}"
            item = json.loads((directory / "OUTPUTS_SEALED.json").read_text())
            if item.get("reference_or_truth_seen") is not False or item.get("outer_TX_seen") is not False:
                raise ValueError("case seal has an exposure violation")
            for filename in ("QUERIES_AND_CONDITIONALS.npz", "ONLINE.json", "QUERY_TRACE.json",
                             "DEFAULT_UNKNOWN_TAIL_INTERVALS.npz", *(m+".npz" for m in METHODS)):
                if not (directory / filename).is_file():
                    raise ValueError("sealed output missing: "+str(directory/filename))
    return seal


def quantiles(value):
    value = np.asarray(value)
    return {name:float(v) for name,v in zip(("minimum","median","p90","maximum"),
                                          np.quantile(value,[0,.5,.9,1]))}


def online_phase(output, old_root, prior, observations, schedule, book):
    output.mkdir()
    write_json(output / "PRIOR_SCHEDULE_BEFORE_QUERIES.json", dict(
        chosen_before_any_response_query=True, prior_id=PRIOR_ID, selected_wall_indices=schedule.selected,
        selected_wall_group_ids=prior.group_ids[schedule.selected], anchor_target_index=schedule.anchor,
        anchor_pose=prior.pose_nodes[schedule.anchor], unique_atom_budget_per_case=schedule.unique_queries,
        rule="one maximum-prior wall per support, then prior-weighted farthest point; stable wallID ties",
        interpolation="same three nearest same-support selected rows, inverse squared metadata distance",
        schedule_independent_of_y_truth_and_response=True))
    write_json(output / "SOURCE_ACCESS.json", dict(prior_and_base_y_only=True,
        data_source="individual paid_rows only; streamed leading TX0/1 bytes",
        FINITE_BANK_opened=False, truth_opened=False, references_opened=False,
        outer_TX_values_loaded=False, root=old_root,
        observed_y_shape=list(observations.shape), prior_id=PRIOR_ID))
    oracle = BaseRowOracle(old_root)
    cases=[]
    for parent in range(12):
        for repeat in range(2):
            before = len(oracle.calls)
            result=online_core(prior, observations[parent,repeat], oracle, schedule)
            count=len(oracle.calls)-before
            save_online(output/f"parent_{parent:02d}_repeat_{repeat:02d}",result,prior)
            cases.append(dict(parent=parent,repeat=repeat,shared_queries=count,status=result["status"],
                sampled_diameter=result["sampled_conditional_TV_diameter"]))
            event(book,"old_bank_online_case",parent=parent,repeat=repeat,
                  logical_atom_queries=count,standalone_per_method_queries=count,
                  shared_by_three_methods=True,actual_new_full_FE=0,actual_new_RHS=0)
            print(f"online parent={parent:02d} repeat={repeat:02d} queries={count} consensus=NOT_CERTIFIED",flush=True)
    write_json(output/"ALL_OUTPUTS_SEALED.json",dict(case_count=24,method_ids=METHODS,
        all_method_outputs_written=True,reference_or_truth_seen=False,outer_TX_seen=False,
        total_shared_unique_atom_queries=len(oracle.calls),cases=cases,label=LABEL))


def score_phase(output, online, old_root, prior, observations, schedule, book):
    verify_seal(online)  # Must precede construction of every teacher/reference and truth read.
    output.mkdir()
    write_json(output/"SOURCE_ACCESS.json",dict(all_online_outputs_verified_before_access=True,
        online_files_modified=False,full_teacher_used_only_after_global_seal=True,
        teacher_source="2950 individual base-port paid_rows, streamed TX0/1 only",
        truth_source=old_root.parent/"offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz",
        truth_keys=["truth_index"],outer_TX_values_loaded=False,
        query_savings_claim=False,label=LABEL))
    oracle=BaseRowOracle(old_root)
    atom_responses=np.array([oracle(atom) for atom in range(2950)])
    response=atom_responses[prior.atom_by_wall_target]
    truth=load_keys(old_root.parent/"offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz",("truth_index",))["truth_index"]
    if truth.shape!=(12,) or truth.dtype.kind not in "iu":
        raise ValueError("twelve original truth atom indices required")
    atom_pose=np.empty((2950,2)); atom_wall=np.empty(2950,dtype=int)
    for wall in range(118):
        atom_pose[prior.atom_by_wall_target[wall]]=prior.pose_nodes
        atom_wall[prior.atom_by_wall_target[wall]]=wall
    static=static_distances(response,prior.central_target)
    np.savez_compressed(output/"OFFLINE_STATIC_DISTANCES.npz",**static)
    event(book,"old_bank_full_teacher_after_seal",logical_atom_queries=2950,
          old_response_cache_setup_queries=2950,actual_new_full_FE=0,actual_new_RHS=0,
          offline_only=True,query_savings_claim=False)
    rows=[]; diag_rows=[]; metric_summaries=[]
    for parent in range(12):
        wall_truth=atom_wall[truth[parent]]
        family="empty_wall" if not np.any(prior.support[wall_truth]) else "nonempty_wall"
        support_family=f"support_count_{int(np.sum(prior.support[wall_truth]))}"
        for repeat in range(2):
            directory=online/f"parent_{parent:02d}_repeat_{repeat:02d}"
            case_meta=json.loads((directory/"ONLINE.json").read_text())
            reference=conditional_reference(prior,response,observations[parent,repeat])
            consensus=load_keys(directory/(METHODS[2]+".npz"),("target_weights",))["target_weights"]
            distance=offline_distances(prior,response,observations[parent,repeat],reference,static,schedule.selected,consensus)
            np.savez_compressed(output/f"OFFLINE_REFERENCE_parent_{parent:02d}_repeat_{repeat:02d}.npz",**reference)
            np.savez_compressed(output/f"OFFLINE_DISTANCES_parent_{parent:02d}_repeat_{repeat:02d}.npz",**distance)
            refq=reference["target_weights"]; refmean=refq@prior.pose_nodes
            refgt=float(np.sum((refmean-atom_pose[truth[parent]])**2))
            refrisk=float(np.sum(refq*np.sum((prior.pose_nodes-refmean)**2,axis=1)))
            intervals=load_keys(directory/"DEFAULT_UNKNOWN_TAIL_INTERVALS.npz",("target_lower","target_upper"))
            for method in METHODS:
                q=load_keys(directory/(method+".npz"),("target_weights",))["target_weights"]
                if q.shape!=(25,) or np.any(q<0) or not np.isclose(q.sum(),1):
                    raise ValueError("sealed proper target posterior required")
                mean=q@prior.pose_nodes
                row=dict(parent=parent,repeat=repeat,family=family,support_family=support_family,
                    method=method,shared_queries=790,target_TV=float(.5*np.abs(q-refq).sum()),
                    position_Bayes_regret=float(np.sum((mean-refmean)**2)),
                    position_GT_squared_error=float(np.sum((mean-atom_pose[truth[parent]])**2)),
                    reference_GT_squared_error=refgt,reference_position_Bayes_risk=refrisk,
                    strict_saved_finite_TV_upper=case_meta["method_metadata"][method]["strict_saved_finite_TV_upper"],
                    target_interval_max_violation=float(max(np.max(intervals["target_lower"]-refq),
                                                            np.max(refq-intervals["target_upper"]))),
                    TV_gate_002_pass=bool(.5*np.abs(q-refq).sum()<=TV_GATE),
                    online_certified=False,actual_new_full_FE=0,actual_new_RHS=0)
                for axis,name in enumerate(("x","y")):
                    row["CDF_"+name+"_sup"]=max(abs(float(q[prior.pose_nodes[:,axis]<=v].sum()
                         -refq[prior.pose_nodes[:,axis]<=v].sum())) for v in np.unique(prior.pose_nodes[:,axis]))
                for name,mask in (("x>=0",prior.pose_nodes[:,0]>=0),("y>=2.1",prior.pose_nodes[:,1]>=2.1)):
                    row["event_"+name+"_signed_error"]=float(q[mask].sum()-refq[mask].sum())
                rows.append(row)
            diag=dict(parent=parent,repeat=repeat,family=family,support_family=support_family,
                **{k:v for k,v in distance.items() if not isinstance(v,np.ndarray)},
                online_sampled_diameter=case_meta["sampled_conditional_TV_diameter"],
                online_sampled_radius=case_meta["sampled_conditional_TV_radius"],
                online_unseen_radius_bound=1.0,consensus_certificate=False,
                reference_mixture_identity_error=reference["mixture_identity_max_error"],
                full_table_diameter_below_002=distance["all_conditional_TV_diameter"]<=TV_GATE)
            diag_rows.append(diag)
            summary=dict(parent=parent,repeat=repeat,scope="OFFLINE_FULL_TABLE_NO_QUERY_SAVINGS_CLAIM")
            for name in ("contrast_residual_aware_Linf","contrast_residual_aware_oscillation","conditional_target_TV","conditional_mean_squared_distance"):
                summary[name]=quantiles(distance[name])
                summary[name+"_sampled_pairs"]=quantiles(distance[name][distance["sampled_pair_mask"]])
            metric_summaries.append(summary)
            event(book,"old_bank_reference_score",parent=parent,repeat=repeat,
                full_reference_logical_atom_evaluations=2950,actual_new_full_FE=0,actual_new_RHS=0,offline_only=True)
    write_csv(output/"TARGET_RESULTS.csv",rows)
    write_csv(output/"CONDITIONAL_DIAGNOSTICS.csv",diag_rows)
    parent_rows=[]
    for parent in range(12):
        for method in METHODS:
            subset=[r for r in rows if r["parent"]==parent and r["method"]==method]
            parent_rows.append(dict(parent=parent,family=subset[0]["family"],support_family=subset[0]["support_family"],method=method,
                repeats=2,target_TV_mean=float(np.mean([r["target_TV"] for r in subset])),
                target_TV_max=max(r["target_TV"] for r in subset),
                position_Bayes_regret_mean=float(np.mean([r["position_Bayes_regret"] for r in subset])),
                position_GT_MSE=float(np.mean([r["position_GT_squared_error"] for r in subset])),
                position_GT_RMSE=float(np.sqrt(np.mean([r["position_GT_squared_error"] for r in subset]))),
                reference_GT_MSE=float(np.mean([r["reference_GT_squared_error"] for r in subset])),
                online_certified=False))
    write_csv(output/"PARENT_RESULTS.csv",parent_rows)
    family_rows=[]
    for grouping in ("all","family","support_family"):
        for family in (["all"] if grouping=="all" else sorted({r[grouping] for r in parent_rows})):
            for method in METHODS:
                subset=[r for r in parent_rows if r["method"]==method and (grouping=="all" or r[grouping]==family)]
                family_rows.append(dict(grouping=grouping,family=family,method=method,parents=len(subset),
                    target_TV_parent_mean=float(np.mean([r["target_TV_mean"] for r in subset])),
                    target_TV_max=max(r["target_TV_max"] for r in subset),
                    position_GT_RMSE=float(np.sqrt(np.mean([r["position_GT_MSE"] for r in subset]))),
                    reference_GT_RMSE=float(np.sqrt(np.mean([r["reference_GT_MSE"] for r in subset]))),
                    position_Bayes_regret_parent_mean=float(np.mean([r["position_Bayes_regret_mean"] for r in subset])),
                    online_certified_parent_count=0))
    write_csv(output/"FAMILY_RESULTS.csv",family_rows)
    write_json(output/"OFFLINE_DISTANCE_SUMMARIES.json",dict(
        scope="FULL_OLD_BANK_TEACHER_AFTER_SEAL_NO_QUERY_SAVINGS_CLAIM",
        complete_response_whitened_real_RMS=quantiles(static["complete_response_whitened_real_RMS"]),
        target_increment_whitened_real_RMS=quantiles(static["target_increment_whitened_real_RMS"]),
        background_projection_norm=quantiles(static["background_difference_projection_into_full_increment_space_norm"]),
        increment_space_empirical_rank=static["full_increment_space_empirical_rank"],cases=metric_summaries))
    summary=dict(label=LABEL,status="LOGIC_REPLAY_COMPLETE_NO_ONLINE_CONSENSUS_CERTIFICATE",
        case_count=24,case_method_rows=len(rows),actual_new_full_FE=0,actual_new_RHS=0,
        shared_online_queries_per_case=790,total_shared_online_queries=24*790,
        standalone_queries_per_method_per_case=790,teacher_unique_setup_queries=2950,
        full_reference_logical_atom_evaluations=24*2950,online_certificates=0,
        sampled_diameter_pass_count=sum(r["online_sampled_diameter"]<=TV_GATE for r in diag_rows),
        all_table_diameter_pass_count=sum(r["full_table_diameter_below_002"] for r in diag_rows),
        sampled_diameter_range=[min(r["online_sampled_diameter"] for r in diag_rows),max(r["online_sampled_diameter"] for r in diag_rows)],
        all_diameter_range=[min(r["all_conditional_TV_diameter"] for r in diag_rows),max(r["all_conditional_TV_diameter"] for r in diag_rows)],
        unqueried_background_posterior_mass_range=[min(r["unqueried_background_posterior_mass"] for r in diag_rows),max(r["unqueried_background_posterior_mass"] for r in diag_rows)],
        cold_PDE_triggered=False,cold_PDE_status="NOT_TRIGGERED_NO_RELIABLE_UNQUERIED_CONDITIONAL_BOUND",
        relative_mass_must_be_retained_where_conditionals_disagree=True,
        old_M4_M5_rerun=False,full_table_diagnostics_count_as_online_savings=False,
        old_bank_setup=dict(full_FE=2950,RHS=11800,CPU_seconds=9288.023039999998),
        floating_point_roundoff_certified=False,method_family_results=family_rows)
    write_json(output/"SUMMARY.json",summary)
    return summary


def write_proof_note(path):
    path.write_text("""# E3 条件共识的范围与证明

此单元只做旧 118 墙 × 25 位置 counting prior 的受控逻辑回放，标签为 `REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE`。目标固有参数已知；不验证未知材料、孔洞、连续目标或冷 PDE 加速。

设条件目标分布为 K_w，背景后验权重为 beta。精确目标分布为 beta K。若 delta=max TV(K_w,K_v)，则 TV(beta K,beta' K) <= delta TV(beta,beta')。证明：令 t=TV(beta,beta')。t=0 时直接成立；否则把 beta-beta' 的正负部分分别除以 t，得到两个概率权重 u,v。其条件混合之差的 TV 不超过 sum_wv u_w v_v TV(K_w,K_v) <= delta，因此原差不超过 t delta。此为一般 Dobrushin 概率工具，不主张原创。所有条件距共同 q 的 TV <= e 时，任意背景混合距 q <= e，无须求 beta。

平方位置损失下，共同动作 a 的条件 regret 为 ||a-E_w X||²。若对所有 w 都 <= e，则任何混合下 regret <= e；只检查采样墙的条件均值，不能把未见墙纳入。后验共识与共同动作均不保证真实目标误差小，必须另报 GT RMSE。

Gaussian 白化模型的精确式为 D_w(x)=<r_w0,Delta_w(x)>-||Delta_w(x)||²/2，其中 r_w0=y-F(w,x0)，Delta_w=F(w,x)-F(w,x0)。若 U 包含全族增量，则只需 P_U r_w0 与完整增量函数；相同增量本身不够。两墙差满足 D_1-D_2=<rbar_x,deltaDelta>-<delta f0,barDelta>。本实现保留全部残差投影项；离线全表增量空间与距离不用于选组或免费界。

背景逐行尺度反例：等 prior，L1=(1,eps)，L2=(eps,1)。将 L2 乘 eta，不改变该行内部 contrast，却改变目标混合。`ALGEBRA_CHECKS.json` 保存数值实例。缺失原子反例固定所有墙 anchor 为 t，并完全查询第一行 (t,t²)；未查询行可为 (t,t²)，也可为 (t,1)。两张表的查询记录完全一致，而目标后验可相差接近 1。因此采样条件一致不能覆盖未查询条件。此是抽象概率表，不冒充 Maxwell 物理反例。

本回放没有可独立获得的未查询物理界；未查询 likelihood 默认 [0,1]。在线的全体条件 TV 半径上界为 1，停止判定始终 NOT_CERTIFIED；即使离线发现接近，也不能追溯升级为免费证书。仅当独立物理界覆盖所有正 prior 背景时，才可重审冷 PDE 入口。
""", encoding="utf-8")


def verify_existing(output,book=None):
    """Check saved artifacts only: no oracle, teacher, truth or new queries."""
    if book is None:
        from r2.accounting import BudgetBook
        book=BudgetBook()
    online=output/"online";scores=output/"offline_scores"
    verify_seal(online)
    schedule=json.loads((online/"PRIOR_SCHEDULE_BEFORE_QUERIES.json").read_text())
    selected=np.asarray(schedule["selected_wall_indices"],dtype=int)
    anchor=schedule["anchor_target_index"]
    maxima=dict(query_reference_log_likelihood_error=0.,selected_conditional_reference_error=0.,
                probability_sum_error=0.,interval_signed_violation=-np.inf,
                reference_mixture_error=0.,TV_bound_signed_violation=-np.inf)
    for parent in range(12):
        for repeat in range(2):
            directory=online/f"parent_{parent:02d}_repeat_{repeat:02d}"
            saved=load_keys(directory/"QUERIES_AND_CONDITIONALS.npz",(
                "query_indices","queried_response_data","queried_log_likelihood",
                "selected_conditionals","selected_walls","pose_nodes"))
            reference=load_keys(scores/f"OFFLINE_REFERENCE_parent_{parent:02d}_repeat_{repeat:02d}.npz",(
                "log_likelihood","conditionals","background_weights","target_weights"))
            trace=json.loads((directory/"QUERY_TRACE.json").read_text())
            meta=json.loads((directory/"ONLINE.json").read_text())
            indices=saved["query_indices"]
            assert len(indices)==790 and len(set(indices.tolist()))==790 and len(trace)==790
            assert saved["queried_response_data"].shape==(790,2,16)
            assert np.all(np.isfinite(saved["queried_response_data"]))
            assert np.array_equal(saved["selected_walls"],selected)
            assert all(t["query_ordinal"]==i+1 and t["atom_index"]==int(indices[i]) for i,t in enumerate(trace))
            assert all(t["target_index"]==anchor for t in trace[:118])
            assert {t["wall_index"] for t in trace[:118]}==set(range(118))
            assert {t["wall_index"] for t in trace[118:]}==set(selected.tolist())
            for wall in selected:
                assert {t["target_index"] for t in trace if t["wall_index"]==wall}==set(range(25))
            known=np.isfinite(saved["queried_log_likelihood"])
            err=float(np.max(np.abs(saved["queried_log_likelihood"][known]-reference["log_likelihood"][known])))
            maxima["query_reference_log_likelihood_error"]=max(maxima["query_reference_log_likelihood_error"],err)
            maxima["selected_conditional_reference_error"]=max(maxima["selected_conditional_reference_error"],
                float(np.max(np.abs(saved["selected_conditionals"]-reference["conditionals"][selected]))))
            refq=reference["target_weights"]
            maxima["reference_mixture_error"]=max(maxima["reference_mixture_error"],
                float(np.max(np.abs(reference["background_weights"]@reference["conditionals"]-refq))))
            intervals=load_keys(directory/"DEFAULT_UNKNOWN_TAIL_INTERVALS.npz",("target_lower","target_upper"))
            maxima["interval_signed_violation"]=max(maxima["interval_signed_violation"],
                float(np.max(intervals["target_lower"]-refq)),float(np.max(refq-intervals["target_upper"])))
            assert meta["online_conditional_consensus_certificate"] is False
            assert meta["all_conditionals_TV_radius_upper_without_physical_bounds"]==1
            for method in METHODS:
                q=load_keys(directory/(method+".npz"),("target_weights",))["target_weights"]
                maxima["probability_sum_error"]=max(maxima["probability_sum_error"],abs(float(q.sum())-1.))
                upper=meta["method_metadata"][method]["strict_saved_finite_TV_upper"]
                maxima["TV_bound_signed_violation"]=max(maxima["TV_bound_signed_violation"],float(.5*np.abs(q-refq).sum()-upper))
    assert maxima["query_reference_log_likelihood_error"]<1e-10
    assert maxima["selected_conditional_reference_error"]<1e-12
    assert maxima["probability_sum_error"]<1e-12
    assert maxima["reference_mixture_error"]<1e-12
    assert maxima["interval_signed_violation"]<1e-12
    assert maxima["TV_bound_signed_violation"]<1e-12
    with (scores/"TARGET_RESULTS.csv").open() as stream:
        target_rows=list(csv.DictReader(stream))
    with (scores/"CONDITIONAL_DIAGNOSTICS.csv").open() as stream:
        diagnostics=list(csv.DictReader(stream))
    assert len(target_rows)==72 and len(diagnostics)==24
    summary=json.loads((scores/"SUMMARY.json").read_text())
    receipt=json.loads((ROOT/"results/jobs/e3-consensus-logic/RECEIPT.json").read_text())
    metric=json.loads((scores/"OFFLINE_DISTANCE_SUMMARIES.json").read_text())
    passes={method:sum(r["TV_gate_002_pass"]=="True" for r in target_rows if r["method"]==method) for method in METHODS}
    identity_max=max(float(r["residual_projection_identity_max_abs_error"]) for r in diagnostics)
    pair_identity_max=max(float(r["two_background_contrast_identity_max_abs_error"]) for r in diagnostics)
    result=dict(status="PASS_SAVED_ARTIFACT_CONSISTENCY_NOT_PHYSICS_CERTIFICATE",cases=24,
        methods=3,unique_queries_per_case=790,new_oracle_queries=0,actual_new_full_FE=0,actual_new_RHS=0,
        maxima=maxima,residual_identity_max_abs_error=identity_max,
        two_background_identity_max_abs_error=pair_identity_max,TV_002_pass_counts=passes,
        no_online_files_modified=True,roundoff_certified=False)
    write_json(output/"VERIFICATION.json",result)
    delegated=ROOT/"research/delegated/consensus"
    all_rows=[r for r in summary["method_family_results"] if r["grouping"]=="all"]
    lines=["# E3 旧 bank 条件目标共识：零新 PDE 回放", "",
        "本轮结论：此旧 bank 的全体条件目标分布明显分歧，条件共识不能支持省略背景后验相对质量。24 个案例均为 `NOT_CERTIFIED`；未触发冷 PDE。标签为 `REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE`。", "",
        "在线只开放有限 prior 元数据、TX0/1 主观测 y 和单原子主端口 API。28 个完整曲线由 prior 元数据预先选定，各 support 分支先取一墙，再按 prior 加权最远点填满；另付全部 118 墙的中心 anchor。三种方法共享同一 790 原子记录。条件共识只用已查询行的条件目标分布均匀混合，不计算背景后验权重；它是有意标为未认证的候选。", "",
        f"24 个案例先全部封存，之后才读取 2950 条旧主端口响应生成 teacher 并读取旧 truth_index 评分。新增 FE/RHS 均为 0；在线共享查询共 18,960，teacher 建立 2,950 原子，24 个参考各 2,950 原子等价评估，共 70,800。历史 bank 建立成本仍为 2,950 FE、11,800 RHS、9,288.02304 CPU 秒。此次主任务 wrapper 为 {receipt['process_cpu_seconds']:.6f} CPU 秒、{receipt['wall_seconds']:.3f} 墙钟秒、峰值 RSS {receipt['polled_peak_RSS_bytes']:,} 字节。未重复 M4/M5，未作新物理加速主张。", "",
        "| 方法 | 平均目标 TV | 最大 TV | TV≤0.02 案例 | GT RMSE | 相对参考均值的平方 regret |", "|---|---:|---:|---:|---:|---:|"]
    for row in all_rows:
        lines.append(f"| {row['method']} | {row['target_TV_parent_mean']:.6f} | {row['target_TV_max']:.6f} | {passes[row['method']]}/24 | {row['position_GT_RMSE']:.6f} | {row['position_Bayes_regret_parent_mean']:.6f} |")
    lines.extend(["",f"全有限参考的 GT RMSE 为 {all_rows[0]['reference_GT_RMSE']:.6f}。先在每个 parent 内平均两次噪声重复，再等权平均 12 个 parent；RMSE 为该 MSE 的平方根。parent 与 family 分解见 `PARENT_RESULTS.csv`、`FAMILY_RESULTS.csv`。empty/nonempty 与 support_count 分组仅在封存后由 truth 产生，未用于查询。", "",
        "采样条件 TV 直径范围为 %.9f–%.9f，全表为 %.9f–%.9f；全体 24 个案例均远大于 0.02。直径出现 1+约 2e-16 属浮点舍入，未声称严格舍入认证。未查询背景的实际后验质量从 %.6g 到 %.9f，因此不能用统一删除未查询质量来救场。" % (*summary['sampled_diameter_range'],*summary['all_diameter_range'],*summary['unqueried_background_posterior_mass_range']), "",
        "完整响应与目标增量的白化实 RMS 配对距离中位数分别为 %.6f 和 %.6f；背景差在全族增量空间内的投影范数中位数为 %.6f。该全表有限增量矩阵在相对奇异值阈值 1e-12 下的经验秩为 %d（实观测维数 64），不能把增量变化较小直接当作残差投影可忽略。四类完整配对距离保存在 OFFLINE_STATIC_DISTANCES.npz 与各 OFFLINE_DISTANCES_*.npz；这些全 teacher 诊断不产生在线查询节省。" % (metric['complete_response_whitened_real_RMS']['median'],metric['target_increment_whitened_real_RMS']['median'],metric['background_projection_norm']['median'],metric['increment_space_empirical_rank']), "",
        "默认未知 likelihood 始终使用 [0,1]，全体条件 TV 半径上界为 1；采样直径不覆盖未查询墙。有限 saved-Gaussian 后验区间与所有已付原子均经过封存后回读校对。完整残差投影恒等式最大绝对差为 %.3g，两背景差分恒等式为 %.3g；一般 Dobrushin 混合证明、逐行尺度反例与未查询原子反例详见 PROOF_AND_SCOPE_ZH.md、ALGEBRA_CHECKS.json。它们是代数与一般概率工具，不是新增 Maxwell 世界或物理证书。" % (identity_max,pair_identity_max), "",
        "保留边界：此结果只否定本冻结表示和旧离散 prior 上的全体条件共识入口。它不否定所有物理分组或局部条件共识路线；后者仍需在完整表之外可获得且覆盖未见背景的物理界。没有该界，不启动冷 PDE。", ""])
    (delegated/"REPORT_ZH.md").write_text("\n".join(lines),encoding="utf-8")
    event(book,"consensus_saved_artifacts_verified",new_oracle_queries=0,actual_new_full_FE=0,actual_new_RHS=0,
          status=result["status"])
    print(json.dumps(plain(result)),flush=True)
    return result


def main(argv=None,book=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-bank-root",type=Path,default=(ROOT/"inputs/old_finite/online" if (ROOT/"inputs/old_finite/online").exists() else OLD/"inputs/old_finite/online"))
    parser.add_argument("--output",type=Path,default=ROOT/"results/e3_logic")
    parser.add_argument("--verify-existing",action="store_true")
    args=parser.parse_args(argv)
    old_root=safe_path(args.old_bank_root); output=args.output.resolve()
    if args.verify_existing:
        if not output.is_relative_to(ROOT):
            raise ValueError("verify only under the R2 campaign root")
        return verify_existing(output,book)
    if not output.is_relative_to(ROOT) or output.exists():
        raise ValueError("fresh output directory under the R2 campaign required")
    if book is None:
        from r2.accounting import BudgetBook
        book=BudgetBook()
    soft,hard=resource.getrlimit(resource.RLIMIT_CPU)
    # The parent phase allows more CPU for possible later work; this zero-PDE
    # logic task itself has the narrower predeclared 900-second hard cap.
    limit=min(900,soft) if soft!=resource.RLIM_INFINITY else 900
    resource.setrlimit(resource.RLIMIT_CPU,(limit,hard))
    started=time.process_time()
    output.mkdir(parents=True)
    write_json(output/"FROZEN_PLAN.json",dict(label=LABEL,method_ids=METHODS,
        unique_queries_per_case=790,all_anchor_rows=118,full_conditional_rows=28,
        planned_case_count=24,CPU_cap_seconds=900,TV_gate=TV_GATE,sigma_real_imag=SIGMA,
        online_inputs=["finite prior metadata","base y TX0/1","per-atom base response API"],
        schedule="prior-only stratified farthest point, 28 full curves; 118 paid anchors",
        consensus_point="uniform mixture of 28 sampled exact conditionals, no background posterior weights",
        strict_unqueried_likelihood_interval=[0,1],unseen_conditionals_TV_radius_upper=1,
        references_and_truth_only_after_all_outputs_sealed=True,outer_TX_access=False,
        no_new_PDE=True,old_M4_M5_not_rerun=True,no_query_savings_claim=True))
    write_json(output/"ALGEBRA_CHECKS.json",algebra_checks())
    prior,observations=read_online_prior(old_root)
    schedule=make_schedule(prior)
    online_phase(output/"online",old_root,prior,observations,schedule,book)
    summary=score_phase(output/"offline_scores",output/"online",old_root,prior,observations,schedule,book)
    write_json(output/"COMPLETION.json",dict(status="COMPLETE_REPLAY",label=LABEL,
        own_process_CPU_seconds=time.process_time()-started,actual_new_full_FE=0,actual_new_RHS=0,
        online_certificates=0,cold_PDE_triggered=False,summary=output/"offline_scores/SUMMARY.json"))
    delegated=ROOT/"research/delegated/consensus";delegated.mkdir(parents=True,exist_ok=True)
    write_proof_note(delegated/"PROOF_AND_SCOPE_ZH.md")
    write_json(delegated/"SUMMARY.json",summary)
    print(json.dumps(plain(summary)),flush=True)
    return summary


if __name__=="__main__":
    main()
