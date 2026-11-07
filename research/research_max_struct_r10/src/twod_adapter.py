"""R10 policy adapter over the actual sibling R8 runner; no new physics.

No execution or public-data loading occurs on import. Call execute(args) only
for a parent-approved frozen task. This module never loads GT or trains a NN.
"""
from copy import deepcopy
from datetime import datetime, timezone
import importlib
import importlib.util
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

STUDY = Path(__file__).resolve().parents[1]
REPO = STUDY.parents[1]
ARMS = ("cie_single", "cie_dual", "global_dual", "a_plus", "a_plus_sourcewise")
ALIASES = {"single": "cie_single", "dual": "cie_dual", "global": "global_dual",
           "A": "a_material", "A+global": "a_plus",
           "A+sourcewise": "a_plus_sourcewise", "direct": "physical_lbfgs"}
RULES = ("fixed_001", "noise_ball_99")
PRIVATE_PARTS = {"data_private", "private", "truth", "gt", "sealed_gt", "sealed_truth"}


def resolve(path, base):
    p = Path(path).expanduser()
    return (p if p.is_absolute() else Path(base) / p).resolve()


def public_path(path, base):
    p = resolve(path, base)
    if PRIVATE_PARTS & {v.lower() for v in p.parts}:
        raise ValueError("Private path is not an R10 public input: " + str(p))
    return p


def stamp(path):
    p = Path(path).resolve()
    s = p.stat()
    return dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def dump(path, value):
    import numpy as np
    def convert(v):
        if isinstance(v, np.ndarray):
            return v.tolist()
        if isinstance(v, np.generic):
            return v.item()
        if isinstance(v, Path):
            return str(v)
        raise TypeError(type(v).__name__)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False, default=convert) + "\n")


def inherited_runner(r8_repo):
    """Load R8 by its real path and reject an already imported other copy."""
    repo = Path(r8_repo).resolve()
    path = repo / "research/research_max_struct_r8/experiments/run_case.py"
    if not path.is_file():
        raise FileNotFoundError(path)
    # R8's own absolute __file__ defines its REPO/public path semantics.
    for name, module in list(sys.modules.items()):
        if name.startswith("research.research_max_struct_r10"):
            continue
        if name.startswith("research.") and getattr(module, "__file__", None):
            if not Path(module.__file__).resolve().is_relative_to(repo):
                raise RuntimeError("A different research checkout was already imported: " + name)
    sys.path.insert(0, str(repo))
    spec = importlib.util.spec_from_file_location("r10_inherited_r8_run_case", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require_iid_raw_noise(public, task, config):
    """Only scalar raw circular-complex sigma is supported; no whitening."""
    import numpy as np
    for scope, mapping in (("public", public), ("task", task), ("config", config),
                           ("common", config.get("common", {}))):
        for key, value in mapping.items():
            name = key.lower()
            if any(word in name for word in ("whiten", "covariance", "precision")):
                # A literal false/null flag can explicitly declare no whitening.
                a = np.asarray(value)
                if value is None or (a.shape == () and a.dtype.kind == "b" and not bool(a)):
                    continue
                raise NotImplementedError(f"Unsupported {scope} noise/whitening field: {key}")
            if any(flag in name for flag in ("heteroscedastic", "nonidentical_noise", "source_normalized", "normalized_observations")) and bool(np.asarray(value).any()):
                raise NotImplementedError(f"Unsupported {scope} nonidentical noise or normalized observations")
            if name in ("noise_is_iid", "iid_noise", "noise_iid") and not bool(np.asarray(value).all()):
                raise NotImplementedError(f"Unsupported {scope} non-iid noise flag")
            if name in ("noise_model", "noise_distribution"):
                if str(np.asarray(value).item()) != "iid_circular_complex":
                    raise NotImplementedError(f"Unsupported {scope} noise model: {value}")
            if name in ("sigma_per_source", "sigma_per_receiver", "noise_sigmas", "noise_variance"):
                raise NotImplementedError(f"Unsupported {scope} non-scalar sigma convention: {key}")
    if "sigma" not in public:
        raise ValueError("noise_ball_99 requires public scalar sigma; noise_radius is not substituted")
    a = np.asarray(public["sigma"])
    if a.size != 1 or np.iscomplexobj(a):
        raise NotImplementedError("Only one real sigma shared by all complex observations is supported")
    sigma = float(a.reshape(-1)[0])
    if not np.isfinite(sigma) or sigma < 0:
        raise ValueError("Public sigma must be finite and nonnegative")
    return sigma


def noise_ball_rank(S, Y, sigma, basis_path=None):
    """Smallest supported prefix meeting one full 99% noise-ball discrepancy.

    E|N[r,s]|²=sigma². One event ||N||F²<=R² controls every orthogonal
    projector, even when the chosen prefix is data dependent. No per-r tests,
    rank-dependent chi-square radii, source normalization, whitening, or GT.
    """
    import numpy as np
    from scipy.linalg import svd
    from scipy.stats import chi2
    S, Y = np.asarray(S, complex), np.asarray(Y, complex)
    if S.ndim != 2 or Y.ndim != 2 or S.shape[0] != Y.shape[0] or Y.size == 0:
        raise ValueError("Expected S=(nr,n), Y=(nr,ns), with nonempty data")
    if not np.isfinite(S).all() or not np.isfinite(Y).all():
        raise ValueError("Nonfinite public operator or observations")
    if not np.isfinite(sigma) or sigma < 0:
        raise ValueError("Invalid circular-complex sigma")
    wall, cpu = time.perf_counter(), time.process_time()
    U, singular, Vh = svd(S, full_matrices=False, check_finite=False)
    tolerance = float(np.finfo(float).eps * max(S.shape) * singular[0]) if singular.size else 0.
    supported = int(np.count_nonzero(singular > tolerance))
    coordinates = U.conj().T @ Y
    energies = np.sum(abs(coordinates) ** 2, axis=1)
    # Positive tail accumulation avoids subtracting nearly equal total energies.
    outside = 0. if U.shape[1] == Y.shape[0] else float(np.linalg.norm(Y - U @ coordinates) ** 2)
    tail = outside + np.r_[np.cumsum(energies[::-1])[::-1], 0.]
    radius2 = float(sigma ** 2 * chi2.ppf(.99, 2 * Y.size) / 2)
    if not np.isfinite(radius2):
        raise ValueError("Nonfinite full noise-ball radius")
    eligible = np.flatnonzero(tail[:supported + 1] <= radius2)
    decision = dict(rule="noise_ball_99", probability=.99,
        noise_convention="independent circular complex; E|n|^2=sigma^2",
        sigma=sigma, complex_observations=int(Y.size), real_noise_dof=2 * int(Y.size),
        full_noise_ball_radius_squared=radius2,
        numerically_supported_rank=supported, numerical_singular_threshold=tolerance,
        supported_rule="singular_i > eps(float64)*max(S.shape)*singular_0; same as R8",
        singular_values=singular.tolist(), discarded_energy_squared_by_rank=tail.tolist(),
        new_receiver_SVD=1, new_physical_LU=0, new_complex_RHS=0,
        selection_wall_seconds=time.perf_counter()-wall, selection_cpu_seconds=time.process_time()-cpu,
        guarantee_scope="one full noise event shared by all projectors; no per-r significance or GT")
    if not eligible.size:
        decision.update(status="NO_SUPPORTED_RANK_MEETS_FULL_NOISE_BALL", rank=None)
    else:
        rank = int(eligible[0])
        if basis_path is not None:
            Vd = Vh[:rank].conj().T
            JD = Vd @ (coordinates[:rank] / singular[:rank, None])
            np.savez_compressed(basis_path, Vd=Vd, JD=JD,
                receiver_left=U, singular_values=singular,
                receiver_rank=np.array(rank), grid_shape=np.array([24,24]))
            decision['receiver_basis_snapshot'] = str(basis_path)
            decision['snapshot_new_SVD'] = 0
        decision.update(status="SELECTED", rank=rank, discarded_energy_squared=float(tail[rank]))
    return decision


def execute(args):
    """Prepare an R10 task, then call unchanged R8 execute; no subprocesses."""
    import numpy as np
    from threadpoolctl import threadpool_limits
    wall, cpu = time.perf_counter(), time.process_time()
    r8_repo = resolve(args.r8_repo or str(REPO), REPO)
    runner = inherited_runner(r8_repo)
    manifest_path = public_path(args.manifest, REPO)
    manifest = json.loads(manifest_path.read_text())
    task = deepcopy(runner.choose_task(manifest, args.task_id))
    original_arm = task["arm"]
    task["arm"] = ALIASES.get(original_arm, original_arm)
    if task["arm"] not in ARMS:
        raise ValueError("Unsupported inherited arm: " + str(original_arm))
    config_ref = args.config or task.get("config") or (manifest.get("config") if isinstance(manifest, dict) else None)
    if not config_ref:
        raise ValueError("An explicit frozen solver config is required")
    config_path = public_path(config_ref, REPO)
    cfg = json.loads(config_path.read_text())
    rule = args.dpcs_rule or task.get("dpcs_rule", "fixed_001")
    if rule not in RULES:
        raise ValueError("Unsupported DPCS rule: " + str(rule))
    # Public R8 file references default to the actual frozen sibling checkout.
    input_base = public_path(task.get("input_base", manifest.get("input_base", str(r8_repo))
                             if isinstance(manifest, dict) else str(r8_repo)), REPO)
    if not isinstance(task.get("public_files"), list) or not task["public_files"]:
        raise ValueError("Task must declare public_files")
    task["public_files"] = [str(public_path(f, input_base)) for f in task["public_files"]]
    init = task.get("initialization", "uniform")
    if init not in ("uniform", "nn", "frozen_nn"):
        raise ValueError("R10 adapter accepts uniform or supplied frozen NN material only")
    if init == "uniform" and task.get("initialization_file"):
        raise ValueError("Do not label a supplied material file as uniform")
    if init != "uniform" and not task.get("initialization_file"):
        raise ValueError("Frozen NN initialization requires initialization_file")
    if task.get("initialization_file"):
        task["initialization_file"] = str(public_path(task["initialization_file"], input_base))
    out = resolve(args.output_dir or task.get("output_dir") or
                  str(STUDY / "RESULTS" / cfg["run_id"] / task["task_id"]), REPO)
    diagnostic=(task.get('experiment_scope')=='EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC' and
                out.is_relative_to(STUDY/'RESULTS/right_projection_diagnostic'))
    conditional=(task.get('experiment_scope')=='FRESH_CONDITIONAL_CONFIRMATION' and
                 out.is_relative_to(STUDY/'RESULTS/conditional_confirmation'))
    if not out.is_relative_to(STUDY / "RESULTS" / "twod") and not diagnostic and not conditional:
        raise ValueError("All R10 outputs must stay inside research/research_max_struct_r10")
    if out.exists():
        raise FileExistsError(out)
    task["output_dir"] = str(out)
    prepared = out.parent / (out.name + "__adapter_inputs")
    prepared.mkdir(parents=True, exist_ok=False)
    provenance = dict(adapter="R10 baseline policy/receipt adapter; no new physical implementation",
        inherited_repository=str(r8_repo), inherited_run_case=stamp(runner.__file__),
        inherited_benchmark_solver=stamp(r8_repo / "research/research_max_struct_r8/src/benchmark_solver.py"),
        inherited_cie_fft=stamp(r8_repo / "research/research_max_struct_r8/src/cie_fft.py"),
        original_manifest=stamp(manifest_path), original_config=stamp(config_path),
        public_files=[stamp(f) for f in task["public_files"]], original_arm=original_arm,
        canonical_arm=task["arm"], public_input_base=str(input_base), hash_checks="NONE")
    if task.get("initialization_file"):
        provenance["frozen_initialization_file"] = stamp(task["initialization_file"])
    preparation = dict(status="PREPARING", truth_access="NONE", scoring="NOT_RUN", provenance=provenance)
    dump(prepared / "ADAPTER_PREPARATION.json", preparation)
    inherited_started = inherited_returned = False
    try:
        effective = deepcopy(cfg)
        common = effective["common"]
        if "dpcs_rank" in common or float(common.get("dpcs_rtol", .001)) != .001:
            raise ValueError("Input config must use the separate inherited fixed .001 control")
        common["dpcs_rtol"] = .001
        if init == "uniform":
            value = float(task.get("initial_material", common["host"]))
            lo, hi = common["material_box"]
            if not np.isfinite(value) or not lo <= value <= hi:
                raise ValueError("Uniform material must be finite and in the common box")
        decision = dict(rule=rule, status="INHERITED_FIXED_RTOL", dpcs_rtol=.001,
                        new_receiver_SVD=0, new_physical_LU=0, new_complex_RHS=0)
        with threadpool_limits(limits=1):
            if task.get("initialization_file"):
                with np.load(task["initialization_file"], allow_pickle=False) as z:
                    if {key.lower() for key in z.files} & {"truth", "gt", "m_true", "m_truth", "mask", "family", "seed", "alpha", "y_clean", "yclean", "ground_truth", "true_material"}:
                        raise ValueError("Private payload in supplied NN initializer")
                    m0 = z["m_init"]
                    lo, hi = common["material_box"]
                    if m0.size != 576 or np.iscomplexobj(m0) or not np.isfinite(m0).all() or np.min(m0) < lo or np.max(m0) > hi:
                        raise ValueError("Frozen NN material must be finite, real, 576-dimensional and in the common box")
            if rule == "noise_ball_99":
                public = runner.load_public(task["public_files"])
                sigma = require_iid_raw_noise(public, task, effective)
                if task["arm"] == "physical_lbfgs":
                    decision.update(status="NOT_APPLICABLE_DIRECT_PHYSICS", sigma=sigma,
                                    explanation="direct arm uses no DPCS; no rank-selection SVD charged")
                else:
                    decision = noise_ball_rank(public["S"], public["Y"], sigma, prepared / "RECEIVER_DPCS.npz")
                    if decision["rank"] is None:
                        dump(prepared / "DPCS_RULE.json", decision)
                        raise ValueError(decision["status"])
                    common.pop("dpcs_rtol")
                    common["dpcs_rank"] = decision["rank"]
            task["dpcs_rule"] = rule
            task["config"] = str(prepared / "EFFECTIVE_CONFIG.json")
            effective["r10_adapter"] = dict(rule=rule, inherited_r8_wrapper=True,
                no_new_physics=True, no_GT=True, source_config=str(config_path))
            dump(prepared / "DPCS_RULE.json", decision)
            dump(prepared / "EFFECTIVE_CONFIG.json", effective)
            dump(prepared / "MANIFEST.json", dict(config=task["config"], tasks=[task]))
            prefix_wall, prefix_cpu = time.perf_counter()-wall, time.process_time()-cpu
            preparation.update(status="PREPARED", dpcs=decision,
                wall_seconds=prefix_wall, cpu_seconds=prefix_cpu)
            dump(prepared / "ADAPTER_PREPARATION.json", preparation)
            if args.prepare_only:
                return preparation, prepared
            inherited_started = True
            # Root may later freeze a separate material-only candidate. This
            # hook leaves arm_config/current PR, metrics, basis, budgets and
            # candidate selection untouched. It is unused by the five baselines.
            custom = task.get('material_readout')
            dispatch_module = importlib.import_module('research.research_max_struct_r8.src.benchmark_solver')
            original_dispatch = dispatch_module.run_auxiliary
            if custom:
                if task['arm'] not in ('a_plus','a_plus_sourcewise'):
                    raise ValueError('Custom material readout requires an explicit A+ current control')
                source = public_path(custom['module_file'], input_base)
                spec = importlib.util.spec_from_file_location('r10_root_material_candidate',source)
                candidate = importlib.util.module_from_spec(spec);spec.loader.exec_module(candidate)
                readout = getattr(candidate,custom['callable'])
                if not callable(readout):raise ValueError('Material readout is not callable')
                prepare_name=custom.get('prepare')
                prepare_fn=getattr(candidate,prepare_name) if prepare_name else None
                if prepare_fn is not None and not callable(prepare_fn):
                    raise ValueError('Candidate prepare entry is not callable')
                provenance['custom_material_readout'] = dict(source=stamp(source),
                    callable=custom['callable'],prepare=prepare_name,label=custom['label'],
                    selection='ROOT_FROZEN_SEPARATE_MANIFEST')
                def material_only_dispatch(model,Y,m0,arm,config):
                    if arm != task['arm']:raise ValueError('Unexpected arm in material-only hook')
                    supplied=dict(config);supplied['material_readout']=readout
                    prepare_wall,prepare_cpu=time.perf_counter(),time.process_time()
                    if prepare_fn is not None:
                        supplied=prepare_fn(model,supplied)
                        if not isinstance(supplied,dict):raise TypeError('prepare(model,config) must return config dict')
                    preparation_cost=dict(CPU_seconds=time.process_time()-prepare_cpu,
                        wall_seconds=time.perf_counter()-prepare_wall)
                    # Called after NUMERICAL_SCHEDULE.json serialization: device
                    # or low-rank matrix objects remain runtime-only config.
                    result=dispatch_module.run_cie_fft(model,Y,m0,supplied)
                    for row in result['history']:
                        if 'iteration' in row:
                            row['state_term_material_semantics']='LEGACY_ONEHOP_DIAGNOSTIC_NOT_CUSTOM_MATERIAL_OBJECTIVE'
                            row['custom_material_objective_field']='material_final_block_loss'
                    result['execution_receipt']['custom_prepare_cost']=preparation_cost
                    result['execution_receipt']['custom_material_label']=custom['label']
                    result['execution_receipt']['material_objective_semantics']='material_final_block_loss is candidate fixed-current objective; state_term_material remains legacy onehop diagnostic'
                    return result
                dispatch_module.run_auxiliary=material_only_dispatch
            try:
                receipt, result_dir = runner.execute(SimpleNamespace(
                    manifest=str(prepared / "MANIFEST.json"), public_list=None,
                    task_id=task["task_id"], config=task["config"], output_dir=str(out)))
            finally:
                dispatch_module.run_auxiliary=original_dispatch
            inherited_returned = True
        # Retain the inherited budget file; export the same paid states with the
        # additional R10 public-preparation CPU charged as a separate prefix.
        if (result_dir / "BUDGET_STATES.npz").is_file():
            with np.load(result_dir / "BUDGET_STATES.npz", allow_pickle=False) as z:
                arrays = {k: z[k].copy() for k in z.files}
            arrays["costs"][:, 0] += prefix_cpu
            columns = [str(v) for v in arrays["columns"]]
            columns[0] = "CPU_seconds_including_R10_preparation"
            arrays["columns"] = np.array(columns)
            np.savez_compressed(result_dir / "R10_BUDGET_STATES.npz", **arrays)
        inherited_git = receipt.pop("git_base", None)
        receipt.update(study="R10", adapter_provenance=provenance,
            inherited_r8_git_base_metadata=inherited_git,
            inherited_r8_wrapper=True, no_new_physical_implementation=True,
            initialization_generation_cost_scope=("Supplied frozen NN file; feature/inference/training creation costs not run or included" if init != "uniform" else "Uniform scalar; no NN creation cost"),
            dpcs_rule=decision, adapter_preparation_cpu_seconds=prefix_cpu,
            adapter_preparation_wall_seconds=prefix_wall,
            r10_total_cpu_seconds=time.process_time()-cpu,
            r10_total_wall_seconds=time.perf_counter()-wall,
            r10_cost_scope="R8 costs retained; R10 totals include preparation and inherited execution. Add neither nested total twice. R10 budget copy adds only preparation CPU; finalization tail is separate.")
        dump(result_dir / "R10_DPCS_RULE.json", decision)
        dump(result_dir / "RECEIPT.json", receipt)
        # No extra propagation: frequencies are the inherited basis's own labels.
        if receipt['status'] != 'SEALED_PUBLIC_PREDICTION':
            return receipt, result_dir
        from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
        with np.load(prepared / "RECEIVER_DPCS.npz", allow_pickle=False) as z:
            vd = z['Vd']
        schedule = json.loads((result_dir / 'NUMERICAL_SCHEDULE.json').read_text())
        frequency_arrays = {}
        basis_metadata = []
        for index, stage in enumerate(schedule['stages']):
            basis = NativeFourierAPCS((24,24),stage['mf'],vd)
            frequency_arrays[f'stage_{index:02d}_frequencies'] = basis.frequencies
            basis_metadata.append(dict(stage=index,mf=stage['mf'],columns=basis.k,
                beta_current=stage['beta_current'],beta_material=stage['beta_material']))
        np.savez_compressed(result_dir / 'BASIS_FREQUENCIES.npz', **frequency_arrays)
        dump(result_dir / 'BASIS_METADATA.json', dict(stages=basis_metadata,
            receiver_rank=decision['rank'], basis='(I-Vd Vd*) IDFT_selected',
            FFT_normalization='ortho', flatten_order='C, (y,x,source)',
            explicit_whitened_roundoff_SVD=False,
            receiver_snapshot=str(prepared / 'RECEIVER_DPCS.npz'),
            stage_joint_semantics='m/J/Q/A are paid actual auxiliary stage endpoints; Q=E+GJ',
            coefficient_transport='exact nested-mask zero padding'))
        receipt.update(r10_total_cpu_seconds=time.process_time()-cpu,
                       r10_total_wall_seconds=time.perf_counter()-wall)
        dump(result_dir / 'RECEIPT.json', receipt)
        return receipt, result_dir
    except Exception as exc:
        preparation.update(status="FAILED_R10_ADAPTER", error=repr(exc),
            wall_seconds=time.perf_counter()-wall, cpu_seconds=time.process_time()-cpu,
            inherited_execution_started=inherited_started,
            inherited_execution_returned=inherited_returned)
        dump(prepared / "ADAPTER_PREPARATION.json", preparation)
        raise
