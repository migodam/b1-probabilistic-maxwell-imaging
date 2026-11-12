"""Diagnostic-only frozen C2 runtime builder and bounded passive observer.

Build sources only; never launches an inversion. Primary predictions must not use
this runtime. CLI: python diagnostic_instrumentation.py --help
"""
from __future__ import annotations

import argparse
import ast
import atexit
from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import sys
import time
import zipfile

CHECKPOINTS = (1, 8, 32, 128, 512, 1024, 1500, 2000, 2500, 3000)
MAX_FULL_STATES = 20
_ACTIVE = None


def _json(value):
    """Convert existing values without computing new numerical diagnostics."""
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "tolist"):
        return _json(value.tolist())
    if isinstance(value, float) and not __import__("math").isfinite(value):
        return {"nonfinite": repr(value)}
    return value


class Recorder:
    """One process, bounded full states; no model/projector/optimizer interface."""

    def __init__(self, out, method, metadata, diagnostic_target_iterations=None):
        self.out = Path(out)
        self.method = method
        self.metadata = _json(metadata)
        targets = [] if diagnostic_target_iterations is None else diagnostic_target_iterations
        if not isinstance(targets, (list, tuple)) or len(targets) > 2 or any(
                isinstance(t, bool) or not isinstance(t, int) or t < 1 for t in targets):
            raise ValueError("diagnostic_target_iterations must contain at most two positive integer indices")
        self.target_iterations = tuple(dict.fromkeys(targets))
        self.reached_target_iterations = set()
        self.native_history = []
        self.dr_history = []
        self.static_arrays = {}
        self.stage_metadata = []
        self.states = {}
        self.slots = {}
        self.sequence = 0
        self.peak_full_states = 0
        self.copy_CPU_seconds = 0.0
        self.export_CPU_seconds = 0.0
        self.export_wall_seconds = 0.0
        self.failed = None
        self.flushed = False
        self.final_context = None
        self.dr_rows = {}
        self.output_file_bytes = {}

    def _copy_state(self, kind, stage, iteration, global_iteration, arrays, scalar):
        started = time.process_time()
        # Replacing the stage's last slot every update ensures a BudgetEnd
        # following this observation still leaves the last complete update.
        keys = [f"{kind}:stage:{stage}:last"]
        first = f"{kind}:stage:{stage}:first"
        if first not in self.slots:
            keys.append(first)
        if global_iteration in CHECKPOINTS:
            keys.append(f"{kind}:fixed:{global_iteration}")
        if global_iteration in self.target_iterations:
            keys.append(f"{kind}:target:{global_iteration}")
        prospective = dict(self.slots)
        sid = self.sequence
        for key in keys:
            prospective[key] = sid
        keep = set(prospective.values())
        if len(keep) > MAX_FULL_STATES:
            raise RuntimeError("Diagnostic bounded-state contract exceeded 20; no eviction invented")
        # Release an unreferenced previous-last snapshot before copying its
        # replacement, so even transient live full states stay within the cap.
        self.states = {key: value for key, value in self.states.items()
                       if key in keep}
        copied = {key: value.copy() if hasattr(value, "copy") else deepcopy(value)
                  for key, value in arrays.items()
                  if value is not None}
        row = dict(kind=kind, stage=stage, iteration=iteration,
                   global_iteration=global_iteration, scalar=_json(scalar),
                   arrays=copied)
        self.states[sid] = row
        self.slots = prospective
        if global_iteration in self.target_iterations:
            self.reached_target_iterations.add(global_iteration)
        self.sequence += 1
        self.peak_full_states = max(self.peak_full_states, len(self.states))
        self.copy_CPU_seconds += time.process_time() - started

    def cie_stage(self, *, stage, JD, Vd, Od, Os, Ls, lo, hi, rho,
                  frequencies, beta_j, beta_m, blocks, sourcewise,
                  receiver_rank, mf):
        started = time.process_time()
        if not self.static_arrays:
            values = dict(JD=JD, Vd=Vd, Od=Od, Os=Os, Ls=Ls,
                          lower=lo, upper=hi)
            if blocks is not None:
                values.update(block_indices=blocks.indices,
                              block_diag_g=blocks.diag_g, block_goff=blocks.goff)
            self.static_arrays = {key: value.copy() for key, value in values.items()}
        self.stage_metadata.append(_json(dict(stage=stage, rho=rho,
            frequencies=frequencies, beta_current=beta_j, beta_material=beta_m,
            sourcewise=sourcewise, receiver_rank=receiver_rank, mf=mf,
            phase="stage_entry_before_first_update")))
        self.copy_CPU_seconds += time.process_time() - started

    def cie_update(self, *, stage, iteration, global_iteration, m_before_readout,
                   m, J, Q, A, grad, direction, pr, alpha, beta_j, beta_m,
                   native_history, operation_counts):
        started = time.process_time()
        self.native_history.append(deepcopy(_json(native_history)))
        self.copy_CPU_seconds += time.process_time() - started
        self._copy_state("cie", stage, iteration, global_iteration,
            dict(m_before_readout=m_before_readout, m=m, J=J, Q=Q, A=A,
                 grad=grad, direction=direction, pr=pr, alpha=alpha),
            dict(beta_current=beta_j, beta_material=beta_m,
                 operation_counts=operation_counts,
                 CPU_seconds=time.process_time(),
                 phase="material_complete_before_native_iteration_callback_and_step_increment"))

    def cie_history(self, row):
        started = time.process_time()
        self.native_history.append(deepcopy(_json(row)))
        self.copy_CPU_seconds += time.process_time() - started

    def dr_update(self, *, iteration, pre_X, pre_Q, cX, cQ, cm, aX, aQ,
                  post_X, post_Q, A_info, projector_cost):
        if self.method != "dr":
            return
        started = time.process_time()
        row = _json(dict(iteration=iteration, constitutive_input_index=iteration-1,
            post_iterator_index=iteration, CPU_seconds=time.process_time(),
            shadow_gap=None, shadow_gap_reason="Not computed by native code at every step",
            A_info=A_info))
        self.dr_rows[iteration] = len(self.dr_history)
        self.dr_history.append(row)
        self.copy_CPU_seconds += time.process_time() - started
        self._copy_state("dr", 0, iteration, iteration,
            dict(pre_X=pre_X, pre_Q=pre_Q, cX=cX, cQ=cQ, cm=cm,
                 reflected_A_X=aX, reflected_A_Q=aQ, post_X=post_X, post_Q=post_Q),
            dict(A_info=A_info, projector_cost=projector_cost,
                 phase="DR_update_complete_before_native_physical_audit"))

    def dr_gap(self, iteration, gap, projector_cost):
        if self.method != "dr":
            return
        started = time.process_time()
        index = self.dr_rows.get(iteration)
        if index is not None:
            self.dr_history[index].update(shadow_gap=_json(gap),
                shadow_gap_reason="Existing native audit-time gap, no extra contraction",
                projector_cost=_json(projector_cost))
        self.copy_CPU_seconds += time.process_time() - started

    def finish(self, *, termination, status, error):
        self.final_context = _json(dict(termination=termination,
                                       runner_status=status, runner_error=error))
        self.flush()

    def flush(self):
        if self.flushed:
            return
        self.flushed = True
        cpu, wall = time.process_time(), time.perf_counter()
        try:
            # NumPy is already imported by the solver. No physics is imported.
            import numpy as np
            self.out.mkdir(parents=True, exist_ok=True)
            for name, rows in (("DIAGNOSTIC_NATIVE_HISTORY.jsonl", self.native_history),
                               ("DIAGNOSTIC_DR_STEPS.jsonl", self.dr_history)):
                with (self.out / name).open("w", encoding="utf-8") as stream:
                    for row in rows:
                        stream.write(json.dumps(row, allow_nan=False) + "\n")
            arrays = {f"static__{key}": value
                      for key, value in self.static_arrays.items()}
            index = []
            for sid, row in sorted(self.states.items()):
                prefix = f"state_{sid:06d}__"
                arrays.update({prefix + key: value for key, value in row["arrays"].items()})
                index.append(dict(snapshot_id=sid, array_prefix=prefix,
                    kind=row["kind"], stage=row["stage"], iteration=row["iteration"],
                    global_iteration=row["global_iteration"], scalar=row["scalar"],
                    fields=list(row["arrays"]),
                    slots=[k for k, value in self.slots.items() if value == sid]))
            np.savez(self.out / "DIAGNOSTIC_STATES.npz", **arrays)
            (self.out / "DIAGNOSTIC_STATE_INDEX.json").write_text(
                json.dumps(dict(states=index, stages=self.stage_metadata), indent=2,
                           allow_nan=False) + "\n", encoding="utf-8")
            self.output_file_bytes = {name: (self.out/name).stat().st_size for name in
                ("DIAGNOSTIC_NATIVE_HISTORY.jsonl", "DIAGNOSTIC_DR_STEPS.jsonl",
                 "DIAGNOSTIC_STATES.npz", "DIAGNOSTIC_STATE_INDEX.json")}
        except Exception as exc:
            self.failed = repr(exc)
            print("DIAGNOSTIC_EXPORT_FAILED: " + self.failed, file=sys.stderr, flush=True)
        self.export_CPU_seconds += time.process_time() - cpu
        self.export_wall_seconds += time.perf_counter() - wall
        receipt = self.receipt()
        try:
            (self.out / "DIAGNOSTIC_RECEIPT.json").write_text(
                json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        except Exception as exc:
            self.failed = repr(exc)
            print("DIAGNOSTIC_RECEIPT_FAILED: " + self.failed, file=sys.stderr, flush=True)

    def receipt(self):
        return dict(role="POST_SCORE_OFFLINE_DIAGNOSTIC",
            primary_selection_changed=False, diagnostic_only=True,
            observer_status="FAILED" if self.failed else "EXPORTED" if self.flushed else "BUFFERED",
            observer_error=self.failed, method=self.method, metadata=self.metadata,
            fixed_global_iterations=list(CHECKPOINTS),
            requested_target_iterations=list(self.target_iterations),
            reached_target_iterations=sorted(self.reached_target_iterations),
            unreached_target_iterations=[t for t in self.target_iterations
                                         if t not in self.reached_target_iterations],
            target_policy="Public primary audit iteration references only; capture if naturally reached, no continuation or extra budget",
            retention="Fixed plus at most two requested checkpoints and first/last completed update per stage; latest replaces prior last",
            max_full_states=MAX_FULL_STATES, retained_full_states=len(self.states),
            peak_full_states=self.peak_full_states,
            completed_CIE_update_rows=sum("iteration" in r for r in self.native_history),
            CIE_native_history_rows=len(self.native_history),
            DR_update_rows=len(self.dr_history), stage_entries=len(self.stage_metadata),
            copy_CPU_seconds=self.copy_CPU_seconds,
            export_CPU_seconds=self.export_CPU_seconds,
            export_wall_seconds=self.export_wall_seconds,
            output_file_bytes=self.output_file_bytes,
            final_context=self.final_context,
            costs_are_nested_in_whole_process_CPU=True,
            timer_scope="Observer copies/conversion and export measured separately; receipt write/final interpreter exit covered only by external whole-child CPU",
            deadline_policy="Original absolute process-CPU checks and 120-second configured cap unchanged; observer/export overhead charged, no compensation",
            array_coverage="No incomplete update, no inner color/root arrays, no recomputed gradient/norm/physics; absent states remain absent",
            atexit_limit="Best effort only: fatal kill or interpreter failure can prevent final export")


def configure(out, method, *, diagnostic_target_iterations=None, **metadata):
    global _ACTIVE
    if _ACTIVE is not None:
        raise RuntimeError("Only one diagnostic recorder per process")
    _ACTIVE = Recorder(out, method, metadata, diagnostic_target_iterations)
    atexit.register(_ACTIVE.flush)
    return _ACTIVE


def _observe(name, *args, **payload):
    if _ACTIVE is not None:
        try:
            getattr(_ACTIVE, name)(*args, **payload)
        except Exception as exc:
            _ACTIVE.failed = repr(exc)
            raise


def capture_cie_stage(**payload):
    _observe("cie_stage", **payload)


def capture_cie_update(**payload):
    _observe("cie_update", **payload)


def capture_cie_history(row):
    _observe("cie_history", row)


def capture_dr_update(**payload):
    _observe("dr_update", **payload)


def capture_dr_gap(iteration, gap, projector_cost):
    _observe("dr_gap", iteration, gap, projector_cost)


def _statements(code):
    nodes = ast.parse(code).body
    for node in nodes:
        node._diagnostic_added = True
    return nodes


def _call_name(node):
    return isinstance(node, ast.Expr) and isinstance(node.value, ast.Call) and \
        isinstance(node.value.func, ast.Attribute) and isinstance(node.value.func.value, ast.Name) and \
        node.value.func.value.id == "history" and node.value.func.attr == "append"


def _native_restoration_check(original, instrumented):
    class Strip(ast.NodeTransformer):
        def visit(self, node):
            if getattr(node, "_diagnostic_added", False):
                return None
            return super().visit(node)
    stripped = Strip().visit(deepcopy(instrumented))
    if ast.dump(original, include_attributes=False) != ast.dump(stripped, include_attributes=False):
        raise ValueError("Diagnostic insertion changed native AST outside marked observer statements")


def transform_cie(source):
    original = ast.parse(source)
    tree = deepcopy(original)
    loop = next((x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "_run_cie"), None)
    if loop is None:
        raise ValueError("Expected inherited R8 def _run_cie context missing")
    hits = dict(update=0, stage_end=0, stage_entry=0)

    class Insert(ast.NodeTransformer):
        def visit_Expr(self, node):
            if not _call_name(node):
                return self.generic_visit(node)
            values = node.value.args
            if len(values) != 1 or not isinstance(values[0], ast.Call) or \
                    not isinstance(values[0].func, ast.Name) or values[0].func.id != "dict":
                raise ValueError("Unexpected R8 history.append payload; inspect source context")
            keys = {k.arg for k in values[0].keywords}
            if "global_iteration" in keys and "PR_plus" in keys:
                hits["update"] += 1
                return [node] + _statements("""
capture_cie_update(stage=stage_index, iteration=it+1, global_iteration=iterations+1,
    m_before_readout=m_before_readout, m=m, J=J, Q=Q, A=A,
    grad=grad, direction=direction, pr=pr, alpha=alpha,
    beta_j=beta_j, beta_m=beta_m, native_history=history[-1],
    operation_counts=dict(cost))
""")
            if "event" in keys and "status" in keys:
                hits["stage_end"] += 1
                return [node] + _statements("capture_cie_history(history[-1])")
            raise ValueError("Unrecognized inherited R8 history.append; refuse broad instrumentation")

        def visit_Assign(self, node):
            # Exact stage J/Q construction already computed by original loop.
            targets = [ast.unparse(t) for t in node.targets]
            value = ast.unparse(node.value)
            if targets == ["Q"] and value == "E + GJD + green.apply(J - JD)":
                hits["stage_entry"] += 1
                return [node] + _statements("""
capture_cie_stage(stage=stage_index, JD=JD, Vd=Vd, Od=Od, Os=Os, Ls=Ls,
    lo=lo, hi=hi, rho=rho, frequencies=frequencies,
    beta_j=beta_j, beta_m=beta_m, blocks=blocks, sourcewise=sourcewise,
    receiver_rank=receiver_rank, mf=stage['mf'])
""")
            return self.generic_visit(node)

    # Restrict transformations to the one known function.
    index = tree.body.index(loop)
    tree.body[index] = Insert().visit(loop)
    expected = dict(update=1, stage_end=1, stage_entry=1)
    if hits != expected:
        raise ValueError(f"R8 observer anchors mismatch: {hits}, expected {expected}; source not supported")
    tree.body[0:0] = _statements(
        "from diagnostic_instrumentation import capture_cie_update, capture_cie_history, capture_cie_stage")
    _native_restoration_check(original, tree)
    ast.fix_missing_locations(tree)
    code = ast.unparse(tree) + "\n"
    compile(code, "diagnostic_cie_fft.py", "exec")
    return code, dict(anchors=hits, native_AST_restoration=True)


def transform_runner(source):
    original = ast.parse(source)
    tree = deepcopy(original)
    main = next((x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "main"), None)
    if main is None:
        raise ValueError("Frozen C2 main() context missing")
    hits = dict(configure=0, DR_loop=0, DR_pre=0, DR_post=0, DR_gap=0, export=0)

    class Insert(ast.NodeTransformer):
        def visit_Expr(self, node):
            text = ast.unparse(node)
            if text == "out.mkdir(parents=True, exist_ok=False)":
                hits["configure"] += 1
                return [node] + _statements("""
from diagnostic_instrumentation import configure, capture_dr_update, capture_dr_gap
_diagnostic_trace = configure(out, algorithm,
    diagnostic_target_iterations=task.get('diagnostic_target_iterations', []),
    root=str(root), task_id=args.task_id, campaign=task['campaign'],
    public_input=str(task['input_file']), scale=cfg.get('scale', 6.0),
    CPU_cap=cap, proposal_fraction=proposal_fraction)
""")
            if text == "dump(out / 'RECEIPT.json', receipt)":
                hits["export"] += 1
                return _statements("""
_diagnostic_trace.finish(termination=stop, status=status, error=error)
receipt['diagnostic_recording'] = _diagnostic_trace.receipt()
receipt['diagnostic_only'] = True
""") + [node]
            return self.generic_visit(node)

        def visit_For(self, node):
            if ast.unparse(node.target) != "k" or ast.unparse(node.iter) != "range(int(cfg['steps']))":
                return self.generic_visit(node)
            hits["DR_loop"] += 1
            # The main loop, not AP/online preparation. Only insert at exact
            # existing PC/update/gap boundaries; no numerical expression changes.
            result = []
            for child in node.body:
                text = ast.unparse(child)
                if isinstance(child, ast.Assign) and \
                        [ast.unparse(t) for t in child.targets] == ["(cx, cq, cm)"] and \
                        ast.unparse(child.value) == "P.material(x, q)":
                    result.extend(_statements("_diagnostic_pre_x, _diagnostic_pre_q = x, q"))
                    hits["DR_pre"] += 1
                result.append(child)
                if isinstance(child, ast.If) and "P.affine_ball(2 * cx - x, 2 * cq - q)" in text and \
                        "x + ax - cx" in text:
                    result.extend(_statements("""
capture_dr_update(iteration=k+1, pre_X=_diagnostic_pre_x, pre_Q=_diagnostic_pre_q,
    cX=cx, cQ=cq, cm=cm, aX=ax, aQ=aq, post_X=x, post_Q=q,
    A_info=P.last_A_info, projector_cost=dict(P.costs))
"""))
                    hits["DR_post"] += 1
                if isinstance(child, ast.If):
                    for pos, sub in list(enumerate(child.body)):
                        if isinstance(sub, ast.Assign) and [ast.unparse(t) for t in sub.targets] == ["gap"] and \
                                ast.unparse(sub.value) == "norm_pair(ax - cx, aq - cq)":
                            child.body[pos+1:pos+1] = _statements("capture_dr_gap(k+1, gap, dict(P.costs))")
                            hits["DR_gap"] += 1
            node.body = result
            return node

    index = tree.body.index(main)
    tree.body[index] = Insert().visit(main)
    expected = dict(configure=1, DR_loop=1, DR_pre=1, DR_post=1, DR_gap=1, export=1)
    if hits != expected:
        raise ValueError(f"Frozen C2 runner anchors mismatch: {hits}, expected {expected}; inspect original contexts")
    _native_restoration_check(original, tree)
    ast.fix_missing_locations(tree)
    code = ast.unparse(tree) + "\n"
    compile(code, "run_experiment.py", "exec")
    return code, dict(anchors=hits, native_AST_restoration=True)


def transform_full_current(source):
    old = "from research.research_max_struct_r8.src import cie_fft as inherited"
    new = "import diagnostic_cie_fft as inherited"
    if source.count(old) != 1:
        raise ValueError("Expected one exact R8 inherited import in frozen full_current.py")
    code = source.replace(old, new, 1)
    compile(code, "full_current.py", "exec")
    return code, dict(private_backend_binding_only=True)


def _read_sources(frozen):
    frozen = Path(frozen).resolve()
    if frozen.is_file() and zipfile.is_zipfile(frozen):
        result = {}
        with zipfile.ZipFile(frozen) as z:
            for name in z.namelist():
                part = Path(name.replace("\\", "/"))
                if part.is_absolute() or ".." in part.parts:
                    raise ValueError("Unsafe runtime archive path")
                if not name.endswith("/") and part.parts and part.parts[0] in ("src", "runtime"):
                    result[part.as_posix()] = z.read(name)
        return result
    if not frozen.is_dir():
        raise ValueError("Frozen C2 root must be a directory or runtime ZIP")
    return {p.relative_to(frozen).as_posix(): p.read_bytes()
            for directory in ("src", "runtime") for p in (frozen/directory).rglob("*")
            if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"}


def build_diagnostic_runtime(frozen_C2_root, new_diag_root, legacy_repo=None, *,
                             validate_only=False):
    """Copy src/runtime and instrument private copies. Never copy data/results.

    Pass the exact deployed legacy repository explicitly. Diagnostic manifests
    should contain absolute public input paths; original relative inputs require
    unchanged public files to be supplied separately under new_diag_root.
    """
    frozen = Path(frozen_C2_root).resolve()
    target = Path(new_diag_root).resolve()
    if target == frozen or target in frozen.parents or frozen in target.parents:
        raise ValueError("New diagnostic root must be separate from frozen source root")
    if not validate_only and target.exists():
        raise FileExistsError("Diagnostic runtime root must not already exist")
    if legacy_repo is None:
        candidates = [frozen.parent, *frozen.parents]
        legacy_repo = next((p for p in candidates if
            (p/"research/research_max_struct_r8/src/cie_fft.py").is_file()), None)
    if legacy_repo is None:
        raise ValueError("--legacy-repo required: C2 archive omits inherited R8/R7 sources")
    legacy = Path(legacy_repo).resolve()
    inherited_path = legacy/"research/research_max_struct_r8/src/cie_fft.py"
    sources = _read_sources(frozen)
    for name in ("src/run_experiment.py", "src/full_current.py", "runtime/field_current.py"):
        if name not in sources:
            raise ValueError("Frozen runtime missing " + name)
    runner, runner_check = transform_runner(sources["src/run_experiment.py"].decode())
    backend, backend_check = transform_cie(inherited_path.read_text(encoding="utf-8"))
    adapter, adapter_check = transform_full_current(sources["src/full_current.py"].decode())
    # Every source is compiled as text only. No solver modules or arrays load.
    for name, raw in sources.items():
        if name.endswith(".py"):
            compile(raw, name, "exec")
    receipt = dict(role="POST_SCORE_OFFLINE_DIAGNOSTIC", numerical_execution=False,
        frozen_C2_source=str(frozen), legacy_repo=str(legacy),
        inherited_loop_source=str(inherited_path), new_root=str(target),
        checkpoints=list(CHECKPOINTS), max_full_states=MAX_FULL_STATES,
        optional_task_checkpoint_field="diagnostic_target_iterations: at most two positive integer public-primary audit references; never GT",
        native_settings_or_selection_changed=False,
        primary_runtime_instrumented=False, diagnostic_overhead_paid=True,
        frozen_source_hash=None, hash_reason="User forbids SHA checks",
        validation=dict(runner=runner_check, backend=backend_check, adapter=adapter_check,
                        all_copied_Python_sources_compile=True),
        qualification="Passive arithmetic-preserving observer changes CPU scheduling/stopping; replay is never primary",
        public_inputs="Not copied. Supply absolute reused public paths or unchanged public data separately.",
        scope="Copies src/runtime only; no RESULTS, queues, GT, manifest or public arrays read/copied")
    if validate_only:
        return receipt
    target.mkdir(parents=True, exist_ok=False)
    for name, raw in sources.items():
        dest = target/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
    for name, value in (("src/run_experiment.py", runner),
                        ("src/full_current.py", adapter),
                        ("src/diagnostic_cie_fft.py", backend)):
        (target/name).write_text(value, encoding="utf-8")
    shutil.copyfile(Path(__file__).resolve(), target/"src/diagnostic_instrumentation.py")
    (target/"DIAGNOSTIC_BUILD.json").write_text(json.dumps(receipt, indent=2)+"\n",
                                                encoding="utf-8")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-c2-root", type=Path, required=True,
                        help="Exact extracted C2 deployment root or frozen runtime ZIP")
    parser.add_argument("--new-diag-root", type=Path, required=True,
                        help="New, separate output runtime directory; never a primary root")
    parser.add_argument("--legacy-repo", type=Path,
                        help="Exact existing inherited repository, not a new implementation")
    parser.add_argument("--validate-only", action="store_true",
                        help="AST restoration/compile checks only; no output runtime written")
    args = parser.parse_args()
    receipt = build_diagnostic_runtime(args.frozen_c2_root, args.new_diag_root,
        args.legacy_repo, validate_only=args.validate_only)
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()

