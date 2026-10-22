"""Truth-free learned-component image actions; no training or physics calls.

These integrate a normalized learned weak/probit density times a declared
strong Gaussian reference. They are neither calibrated posterior means nor
teachers. The defensive prior mixture is deliberately not sampled here.
"""
from pathlib import Path
from contextlib import contextmanager
import json

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import ndtr

from p1_material.prior import saturation
from .learned_proposal import LearnedFullProposal

FAMILIES = ("mlp_diag", "mlp_full", "mlp_mdn", "mlp_flow", "cnn_full")
SEEDS = (2026102270, 2026102271, 2026102272)


def validation_selection(out):
    """Require all fifteen NN runs and fitted simple controls; validation only."""
    out = Path(out)
    candidates, reports = [], {}
    for family in FAMILIES:
        for seed in SEEDS:
            closure = out / f"D/jobs/train_{family}_s{seed}/CLOSED.json"
            if json.loads(closure.read_text()).get("status") != "COMPLETE":
                raise RuntimeError("all fifteen training jobs must finish before locked selection")
            folder = out / f"D/models/{family}/seed{seed}"
            report = json.loads((folder / "TRAINING_REPORT.json").read_text())
            metadata = report["metadata"]
            nll = metadata["best_validation_raw_target_nll"]
            checkpoint = (folder / "BEST_DENSITY.pt").resolve()
            if (not metadata.get("best_checkpoint_available") or nll is None or
                    not np.isfinite(nll) or not checkpoint.is_file() or
                    metadata.get("seed") != seed or metadata.get("supervision") != "JOINT_SIMULATION"):
                raise RuntimeError("completed training report/checkpoint contract failed")
            logical_checkpoint = Path(f"D/models/{family}/seed{seed}/BEST_DENSITY.pt")
            original_path = Path(report["checkpoint_path"])
            if tuple(original_path.parts[-len(logical_checkpoint.parts):]) != logical_checkpoint.parts:
                raise RuntimeError("checkpoint logical suffix differs from completed report")
            reports[f"{family}/seed{seed}"] = report
            candidates.append(dict(model_type=family, backend="torch", seed=seed, checkpoint_path=logical_checkpoint.as_posix(),
                                   best_validation_raw_target_nll=float(nll)))
    simple_closure = out / "D/jobs/fit_simple_gaussian_controls/CLOSED.json"
    if json.loads(simple_closure.read_text()).get("status") != "COMPLETE":
        raise RuntimeError("simple Gaussian control fit must finish before locked selection")
    simple = json.loads((out / "D/simple_density/RESULTS.json").read_text())
    if (simple.get("status") != "FITTED_NORMALIZED_JOINT_SUPERVISION_CONTROLS" or
            simple.get("train_parents") != 512 or simple.get("validation_parents") != 128 or
            simple.get("locked_read") is not False):
        raise RuntimeError("simple control result provenance invalid")
    rows = simple["rows"]
    if not rows or any(not np.isfinite(r["validation_raw_NLL"]) for r in rows):
        raise RuntimeError("complete finite simple fit results required")
    if simple["selected"] != min(rows, key=lambda r: r["validation_raw_NLL"]):
        raise RuntimeError("declared simple selection differs from validation minimum")
    unconditional = [r for r in rows if r["dimension"] == 0]
    if simple["unconditional_selected"] != min(unconditional, key=lambda r: r["validation_raw_NLL"]):
        raise RuntimeError("declared unconditional selection differs from validation minimum")
    for label, simple_row in (("simple_conditional_gaussian", simple["selected"]),
                              ("simple_unconditional_gaussian", simple["unconditional_selected"])):
        logical_checkpoint = Path(simple_row["checkpoint"])
        if logical_checkpoint.is_absolute() or ".." in logical_checkpoint.parts or logical_checkpoint.parts[:3] != ("D", "simple_density", "fits"):
            raise RuntimeError("simple checkpoint must have portable D/simple_density/fits path")
        checkpoint = (out / logical_checkpoint).resolve()
        if not checkpoint.is_file():
            raise RuntimeError("selected simple control checkpoint missing")
        candidates.append(dict(model_type=label, backend="numpy", seed=None,
                               checkpoint_path=logical_checkpoint.as_posix(), fit_id=simple_row["id"],
                               best_validation_raw_target_nll=float(simple_row["validation_raw_NLL"])))
    selected = min(candidates, key=lambda r: (r["best_validation_raw_target_nll"], r["model_type"], -1 if r["seed"] is None else r["seed"]))
    return dict(selected=selected, candidates=candidates,
                neural_candidate_reports=reports, simple_candidate_report=simple,
                rule="VALIDATION_NLL_ALL15_AND_SIMPLE_ONLY; ties lexical family then seed; no held-out reads")


def load_selected_density(selected, *, out=None, rng=None):
    """Restore either frozen CPU neural wrapper or pure NumPy simple density."""
    if out is None:
        out = Path(__file__).resolve().parents[2] / "research/p3l_local_posterior_imaging"
    root = Path(out).resolve()
    logical_checkpoint = Path(selected["checkpoint_path"])
    if logical_checkpoint.is_absolute() or ".." in logical_checkpoint.parts:
        raise ValueError("selected checkpoint must be portable relative to campaign OUT")
    checkpoint = (root / logical_checkpoint).resolve()
    if not checkpoint.is_relative_to(root):
        raise ValueError("checkpoint escaped campaign root")
    if selected["backend"] == "torch":
        loaded = load_trained_density(checkpoint)
        if (loaded["model_type"] != selected["model_type"] or
                loaded["metadata"].get("seed") != selected["seed"] or
                loaded["best_validation_raw_target_nll"] != selected["best_validation_raw_target_nll"]):
            raise ValueError("loaded neural family/seed/validation metadata differs from frozen selection")
        return loaded
    if selected["backend"] != "numpy" or selected["model_type"] not in (
            "simple_conditional_gaussian", "simple_unconditional_gaussian"):
        raise ValueError("unsupported selected density backend/family")
    from .simple_density import SimpleConditionalGaussian
    results = json.loads((root / "D/simple_density/RESULTS.json").read_text())
    matching = [row for row in results["rows"] if row["id"] == selected["fit_id"]]
    if (len(matching) != 1 or matching[0]["checkpoint"] != logical_checkpoint.as_posix() or
            matching[0]["validation_raw_NLL"] != selected["best_validation_raw_target_nll"] or
            selected["seed"] is not None):
        raise ValueError("loaded simple fit identity/path/NLL differs from frozen selection")
    if (selected["model_type"] == "simple_unconditional_gaussian") != (matching[0]["dimension"] == 0):
        raise ValueError("simple conditional/unconditional family differs from fitted dimension")
    model = SimpleConditionalGaussian.load(checkpoint, rng=rng)
    dimension = matching[0]["dimension"]
    if (model.arrays["projection"].shape != (1154, dimension) or
            model.arrays["weights"].shape != (dimension, 10)):
        raise ValueError("simple checkpoint PCA/label dimensions differ from fit metadata")
    return dict(model=model, model_type=selected["model_type"], metadata=dict(seed=None),
                best_validation_raw_target_nll=selected["best_validation_raw_target_nll"])


def load_trained_density(checkpoint_path):
    """Inference-only restoration of the stored wrapper, without trainer import."""
    import torch
    from torch import nn
    from .neural import build_model

    class ScaledInferenceDensity(nn.Module):
        def __init__(self, base, scaling):
            super().__init__()
            self.base = base
            for name in ("scalar_mean", "scalar_scale", "label_mean", "label_scale"):
                self.register_buffer(name, torch.as_tensor(scaling[name], dtype=torch.float32).clone())

        def sample(self, count, scalar, image):
            scalar = (scalar - self.scalar_mean) / self.scalar_scale
            return self.base.sample(count, scalar, image) * self.label_scale + self.label_mean

        def log_prob(self, values, scalar, image):
            scalar = (scalar - self.scalar_mean) / self.scalar_scale
            labels = (values - self.label_mean) / self.label_scale
            return self.base.log_prob(labels, scalar, image) - torch.log(self.label_scale).sum()

    payload = torch.load(Path(checkpoint_path), map_location="cpu", weights_only=True)
    if payload.get("format") != "P3L_SUPPLIED_DATA_DENSITY_V1":
        raise ValueError("unsupported trained density format")
    model = ScaledInferenceDensity(build_model(payload["model_type"], **payload["model_options"]), payload["train_scaling"])
    model.load_state_dict(payload["state_dict"], strict=True)
    model.eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return dict(model=model, metadata=payload["metadata"], model_type=payload["model_type"],
                best_validation_raw_target_nll=float(payload["best_validation_raw_target_nll"]))


def context_features(prior, q, prediction, J, observation, sigma):
    """Exactly the training-context scalar and image API; no field solves."""
    q, J = np.asarray(q), np.asarray(J)
    def white(value):
        return np.concatenate((value.real.T, value.imag.T), axis=1).ravel() * np.sqrt(2.) / sigma
    singular = np.linalg.svd(J[:, :64], compute_uv=False)
    scalar = np.r_[white(observation), white(prediction-observation), q, np.log1p(singular*singular)]
    chi = prior.decode(q)[0]
    Dg = prior.derivative(q)
    image = np.stack((chi.real, chi.imag, prior.A @ q[:64],
                      np.sqrt(np.sum(np.abs(Dg[:, :64])**2, axis=1)))).reshape(4, 24, 24)
    if scalar.shape != (1154,) or image.shape != (4, 24, 24):
        raise ValueError("training feature shape mismatch")
    return scalar, image


def make_proposal(model, scalar, image, frame, q0, J, residual):
    proposal = LearnedFullProposal(model, scalar, image, frame["Vw"], frame["Vs"], q0, J, residual)
    if (proposal.d, proposal.k, proposal.s, proposal.p) != (64, 8, 56, 2):
        raise ValueError("frozen full56/weak8/probit2 chart required")
    return proposal


def decode_batch(prior, x):
    """Original decoder formula in ≤512-draw batches, exact probit properties."""
    x = np.asarray(x, float)
    if x.ndim != 2 or x.shape[1] != 66 or not 1 <= len(x) <= 512 or not np.all(np.isfinite(x)):
        raise ValueError("finite (1..512,66) transformed draws required")
    u = ndtr(x[:, 64:])
    real = prior.real_interval[0] + np.diff(prior.real_interval)[0] * u[:, 0]
    imag = prior.imag_interval[0] + np.diff(prior.imag_interval)[0] * u[:, 1]
    field = prior.config["field_mean"] + x[:, :64] @ prior.A.T
    return prior.host + (real+1j*imag)[:, None] * prior.mask[None, :] * saturation(field/prior.temperature)


@contextmanager
def _sampling_stream(proposal, seed):
    if proposal._torch is not None:
        torch = proposal._torch
        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(int(seed))
            yield
    else:
        if not hasattr(proposal.network, "rng"):
            raise ValueError("NumPy density must expose its caller-owned rng")
        previous = proposal.network.rng
        proposal.network.rng = np.random.default_rng(np.random.SeedSequence([int(seed), 0]))
        try:
            yield
        finally:
            proposal.network.rng = previous


def learned_image_actions(proposal, prior, *, count=2048, batch_size=512, seed):
    """IID model-reference MC, using common outer draws for both controls.

    Raw outer values and both complete transformed draw arrays suffice to
    replay material integration. MCSE measures finite IID model draws only;
    it excludes approximation bias, calibration and physical uncertainty.
    """
    if count != 2048 or not isinstance(batch_size, int) or not 1 <= batch_size <= 512:
        raise ValueError("frozen 2048 draws and integer batch size ≤512 required")
    rng = np.random.default_rng(seed)
    accumulators = {name: dict(sum=np.zeros(576, complex), sum_abs2=np.zeros(576))
                    for name in ("full_learned_gaussian", "conditional_mode_plugin")}
    outers, full_draws, plugin_draws = [], [], []
    with _sampling_stream(proposal, seed):
        for left in range(0, count, batch_size):
            n = min(batch_size, count-left)
            # Direct normalized learned-component draw; never proposal.sample(),
            # which intentionally includes the defensive prior mixture.
            outer = proposal._network_call("sample", n)[0]
            if outer.shape != (n, 10) or not np.all(np.isfinite(outer)):
                raise ValueError("learned outer draw shape/finite contract failed")
            a, eta = outer[:, :8], outer[:, 8:]
            b_mode = proposal.strong_mean(a, ndtr(eta))
            noise = solve_triangular(proposal.chol_H.T, rng.normal(size=(56, n)), lower=False).T
            x_full = proposal.from_chart(a, b_mode+noise, eta)
            x_plugin = proposal.from_chart(a, b_mode, eta)
            for name, x in (("full_learned_gaussian", x_full), ("conditional_mode_plugin", x_plugin)):
                material = decode_batch(prior, x)
                accumulators[name]["sum"] += material.sum(axis=0)
                accumulators[name]["sum_abs2"] += np.sum(np.abs(material)**2, axis=0)
            outers.append(outer); full_draws.append(x_full); plugin_draws.append(x_plugin)
    arrays = dict(outer_a8_eta2=np.concatenate(outers), full_x_draws=np.concatenate(full_draws),
                  conditional_mode_x_draws=np.concatenate(plugin_draws))
    diagnostics = {}
    for name, moments in accumulators.items():
        mean = moments["sum"] / count
        raw_variance = (moments["sum_abs2"] - count*np.abs(mean)**2) / (count-1)
        sample_variance = np.maximum(0., raw_variance)
        mcse = np.sqrt(sample_variance/count)
        arrays.update({name+"_material_mean": mean, name+"_material_sum": moments["sum"],
                       name+"_material_sum_abs2": moments["sum_abs2"],
                       name+"_iid_model_mcse_modulus": mcse})
        diagnostics[name] = dict(draws=count, mean_pixel_mcse=float(np.mean(mcse)),
                                 maximum_pixel_mcse=float(np.max(mcse)),
                                 min_raw_sample_variance=float(raw_variance.min()))
    diagnostics.update(seed=int(seed), batch_size=batch_size, defensive_mixture_sampled=False,
        interpretation="LEARNED_OUTER_TIMES_FULL56_GAUSSIAN_REFERENCE; MCSE_ONLY_FINITE_IID_MODEL_DRAWS; NOT_CALIBRATED_POSTERIOR_OR_TEACHER")
    return arrays, diagnostics
