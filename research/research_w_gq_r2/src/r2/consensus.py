"""Access-only finite conditional-consensus logic; no files, bank, truth or PDE.

Sampled conditionals are diagnostics.  A nonzero-prior unqueried conditional
has TV radius upper bound one unless the caller has an independent physical
bound.  This replay deliberately has no such bound or safe early-stop claim.
"""
from dataclasses import dataclass
import numpy as np

from gq.contrast import (STRICT, MarginalEnvelope, logsumexp,
                         normalized_target_intervals, target_posterior)

LABEL = "REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE"
METHODS = ("ORDINARY_RESPONSE_INTERPOLATION", "ANCHOR_WEIGHTED_CONTRAST",
           "SAMPLED_CONDITIONAL_CONSENSUS")
SIGMA = 0.01
TV_GATE = 0.02


def log_array(value):
    with np.errstate(divide="ignore"):
        return np.log(np.asarray(value, dtype=float))


def conditional_prior(prior):
    joint = prior.probability.joint_mass
    wall_mass = joint.sum(axis=1)
    if np.any(wall_mass <= 0):
        raise ValueError("remove undeclared zero-prior wall rows explicitly")
    return wall_mass, joint / wall_mass[:, None]


def potentials(response, y):
    difference = (np.asarray(response) - np.asarray(y)) / SIGMA
    return -0.5 * np.sum(difference.real ** 2 + difference.imag ** 2, axis=(-2, -1))


def real_white(response):
    response = np.asarray(response)
    return np.concatenate((response.real, response.imag), axis=-1).reshape(
        response.shape[:-2] + (-1,)) / SIGMA


def residual_contrast(response, y, anchor):
    """D=<y-f0, f-f0>_R - ||f-f0||_R^2/2, including 1/sigma.

    Complex data are realified: <a,b>_R = Re(a* b), and both real and
    imaginary components have sigma .01.  No background scale is dropped in
    marginal inference; scale disappears only inside each normalized row.
    """
    f = real_white(response)
    yw = real_white(y)
    delta = f - f[:, [anchor]]
    residual0 = yw - f[:, anchor]
    return np.einsum("wm,wxm->wx", residual0, delta) - 0.5 * np.sum(delta ** 2, axis=-1)


def pair_tv(rows):
    rows = np.asarray(rows, dtype=float)
    left, right = np.triu_indices(len(rows), 1)
    distances = 0.5 * np.sum(np.abs(rows[left] - rows[right]), axis=1)
    return (float(distances.max()) if distances.size else 0.0), distances


@dataclass(frozen=True)
class Schedule:
    selected: np.ndarray
    interpolation_weights: np.ndarray
    anchor: int
    unique_queries: int


def make_schedule(prior, count=28):
    """Prior-only stratified farthest-point design, frozen before any query.

    First cover each support branch with its largest prior-mass wall (ID ties).
    Then maximize nearest descriptor distance times sqrt(prior wall mass).
    Response and contrast interpolants use the same three nearest same-support
    selected walls and inverse squared descriptor distance.  Exact rows and
    every paid anchor are retained.
    """
    nw, nx = prior.atom_by_wall_target.shape
    if count != 28 or (nw, nx) != (118, 25):
        raise ValueError("frozen old-bank schedule is 28 of 118 wall curves and 25 targets")
    wall_mass, _ = conditional_prior(prior)
    descriptors = np.asarray(prior.descriptor)
    selected = []
    for branch in np.unique(prior.support_id):
        candidates = np.flatnonzero(prior.support_id == branch)
        selected.append(min(candidates, key=lambda w: (-float(wall_mass[w]), int(prior.group_ids[w]))))
    if len(selected) > count:
        raise ValueError("more support branches than the fixed query budget permits")
    while len(selected) < count:
        def key(wall):
            distance = min(float(np.linalg.norm(descriptors[wall] - descriptors[s])) for s in selected)
            return -distance * np.sqrt(wall_mass[wall]), int(prior.group_ids[wall])
        selected.append(min((w for w in range(nw) if w not in selected), key=key))
    selected = np.asarray(selected, dtype=int)
    weights = np.zeros((nw, count))
    for wall in range(nw):
        if wall in selected:
            weights[wall, np.flatnonzero(selected == wall)[0]] = 1.0
            continue
        candidates = [j for j, s in enumerate(selected) if prior.support_id[s] == prior.support_id[wall]]
        if not candidates:
            raise ValueError("predeclared schedule must cover every support branch")
        neighbors = sorted(candidates, key=lambda j: (
            float(np.linalg.norm(descriptors[wall] - descriptors[selected[j]])),
            int(prior.group_ids[selected[j]])))[:3]
        distances = np.array([np.linalg.norm(descriptors[wall] - descriptors[selected[j]]) for j in neighbors])
        if np.any(distances == 0):
            weights[wall, neighbors[int(np.argmin(distances))]] = 1.0
        else:
            inverse = distances ** -2
            weights[wall, neighbors] = inverse / inverse.sum()
    return Schedule(selected, weights, prior.central_target, nw + count * (nx - 1))


def default_intervals(prior, known_ll):
    """Saved finite Gaussian likelihood [0,1] at unqueried cells only.

    This is the default tail envelope, not a rerun of the old M4 search.  It
    excludes solver/model/roundoff error and cannot supply a Maxwell bound.
    """
    known = np.isfinite(known_ll)
    log_conditional = log_array(prior.probability.wall_given_target)
    lower = logsumexp(log_conditional + np.where(known, known_ll, -np.inf), axis=0)
    upper = logsumexp(log_conditional + np.where(known, known_ll, 0.0), axis=0)
    envelope = MarginalEnvelope(lower, upper, STRICT,
        "default finite saved-Gaussian likelihood [0,1] for each unqueried atom")
    probabilities = normalized_target_intervals(prior.probability.target_mass, envelope)
    atom_intervals = [probabilities.event(mask) for mask in np.eye(len(lower), dtype=bool)]
    events = {}
    for name, selected in (("x>=0", prior.pose_nodes[:, 0] >= 0),
                           ("y>=2.1", prior.pose_nodes[:, 1] >= 2.1)):
        result = probabilities.event(selected)
        events[name] = {"lower": result.lower, "upper": result.upper}
    return dict(log_marginal_lower=lower, log_marginal_upper=upper,
                target_lower=np.array([p.lower for p in atom_intervals]),
                target_upper=np.array([p.upper for p in atom_intervals]), events=events)


def online_core(prior, y, oracle, schedule):
    """The entire online algorithm sees only proper prior, y and query(atom)."""
    y = np.asarray(y, dtype=complex)
    if y.shape != (2, 16) or np.any(~np.isfinite(y)):
        raise ValueError("only finite TX0/1 by 16RX observations are allowed")
    nw, nx = prior.atom_by_wall_target.shape
    response = np.full((nw, nx, 2, 16), np.nan + 1j * np.nan)
    ll = np.full((nw, nx), np.nan)
    indices, paid, traces = [], [], []
    def query(wall, target, phase):
        if np.isfinite(ll[wall, target]):
            return
        atom = int(prior.atom_by_wall_target[wall, target])
        value = np.asarray(oracle(atom), dtype=complex)
        if value.shape != (2, 16) or np.any(~np.isfinite(value)):
            raise ValueError("one finite base-port response per paid atom required")
        response[wall, target] = value
        ll[wall, target] = potentials(value, y)
        indices.append(atom)
        paid.append(value.copy())
        traces.append(dict(query_ordinal=len(indices), atom_index=atom, wall_index=wall,
                           wall_group_id=int(prior.group_ids[wall]), target_index=target,
                           phase=phase, log_likelihood=float(ll[wall, target])))
    for wall in range(nw):
        query(wall, schedule.anchor, "PREDECLARED_ALL_ANCHORS")
    for wall in schedule.selected:
        for target in range(nx):
            query(int(wall), target, "PREDECLARED_PRIOR_ONLY_FULL_CURVE")
    if len(indices) != schedule.unique_queries or len(set(indices)) != len(indices):
        raise ValueError("the frozen shared unique-atom budget differs")
    selected = schedule.selected
    selected_response = response[selected]
    selected_ll = ll[selected]
    contrast = selected_ll - selected_ll[:, [schedule.anchor]]
    contrast_identity = residual_contrast(selected_response, y, schedule.anchor)
    identity_error = float(np.max(np.abs(contrast - contrast_identity)))
    ordinary_response = np.einsum("ws,sxtr->wxtr", schedule.interpolation_weights, selected_response)
    ordinary_response[:, schedule.anchor] = response[:, schedule.anchor]
    ordinary_ll = potentials(ordinary_response, y)
    interpolated_d = schedule.interpolation_weights @ contrast
    contrast_ll = ll[:, [schedule.anchor]] + interpolated_d
    log_wall_given_target = log_array(prior.probability.wall_given_target)
    ordinary_marginal = logsumexp(log_wall_given_target + ordinary_ll, axis=0)
    contrast_marginal = logsumexp(log_wall_given_target + contrast_ll, axis=0)
    ordinary_q, _ = target_posterior(prior.probability.target_mass, ordinary_marginal)
    contrast_q, _ = target_posterior(prior.probability.target_mass, contrast_marginal)
    _, target_given_wall = conditional_prior(prior)
    row_logp = log_array(target_given_wall[selected]) + selected_ll
    row_q = np.exp(row_logp - logsumexp(row_logp, axis=1)[:, None])
    row_q /= row_q.sum(axis=1)[:, None]
    consensus_q = row_q.mean(axis=0)
    consensus_q /= consensus_q.sum()
    diameter, _ = pair_tv(row_q)
    radius = float(np.max(0.5 * np.sum(np.abs(row_q - consensus_q), axis=1)))
    means = row_q @ prior.pose_nodes
    consensus_mean = consensus_q @ prior.pose_nodes
    sampled_regret = float(np.max(np.sum((means - consensus_mean) ** 2, axis=1)))
    roi_diameter_squared = float(np.sum(np.ptp(prior.pose_nodes, axis=0) ** 2))
    intervals = default_intervals(prior, ll)
    outputs = {}
    for method, q in zip(METHODS, (ordinary_q, contrast_q, consensus_q)):
        coordinate_bound = 0.5 * np.sum(np.maximum(np.abs(q - intervals["target_lower"]),
                                                  np.abs(intervals["target_upper"] - q)))
        outputs[method] = dict(target_weights=q, mean=q @ prior.pose_nodes,
            strict_saved_finite_TV_upper=float(min(1.0, coordinate_bound)))
    return dict(label=LABEL, status="NOT_CERTIFIED_UNQUERIED_CONDITIONALS",
                query_indices=np.array(indices, dtype=int), queried_response_data=np.array(paid),
                queried_log_likelihood=ll, query_trace=traces, methods=outputs,
                selected_conditionals=row_q, selected_conditional_means=means,
                selected_contrasts=contrast, selected_walls=selected.copy(),
                interpolation_weights=schedule.interpolation_weights.copy(), intervals=intervals,
                sampled_conditional_TV_diameter=diameter, sampled_conditional_TV_radius=radius,
                sampled_decision_regret_upper=sampled_regret,
                all_conditionals_TV_radius_upper_without_physical_bounds=1.0,
                all_conditionals_decision_regret_upper_without_physical_bounds=roi_diameter_squared,
                sampled_diameter_passes_002=diameter <= TV_GATE,
                online_conditional_consensus_certificate=False,
                online_decision_consensus_certificate=False,
                background_posterior_weights_computed_by_consensus=False,
                unqueried_wall_count=nw-len(selected), unqueried_atom_count=nw*nx-len(indices),
                residual_projection_identity_max_abs_error=identity_error,
                actual_new_full_FE=0, actual_new_RHS=0, safe_physical_query_savings=0,
                full_enumeration_reference_seen=False, floating_point_roundoff_certified=False)


def conditional_reference(prior, response, y):
    """Full old-bank OFFLINE diagnostic; never call before the global seal."""
    ll = potentials(response, y)
    wall_prior, target_given_wall = conditional_prior(prior)
    joint_logs = log_array(target_given_wall) + ll
    row_evidence = logsumexp(joint_logs, axis=1)
    row_q = np.exp(joint_logs - row_evidence[:, None])
    row_q /= row_q.sum(axis=1)[:, None]
    beta_log = log_array(wall_prior) + row_evidence
    beta = np.exp(beta_log - logsumexp(beta_log))
    beta /= beta.sum()
    marginal = logsumexp(log_array(prior.probability.wall_given_target) + ll, axis=0)
    q, logq = target_posterior(prior.probability.target_mass, marginal)
    return dict(log_likelihood=ll, conditionals=row_q, background_weights=beta,
                target_weights=q, log_target_weights=logq, log_marginal=marginal,
                mixture_identity_max_error=float(np.max(np.abs(beta @ row_q - q))))


def static_distances(response, anchor):
    """Full teacher response/increment geometry, charged and labeled offline."""
    f = real_white(response)
    delta = f - f[:, [anchor]]
    left, right = np.triu_indices(len(f), 1)
    response_rms, increment_rms = np.empty(len(left)), np.empty(len(left))
    for start in range(0, len(left), 128):
        part = slice(start, start + 128)
        response_rms[part] = np.sqrt(np.mean((f[left[part]] - f[right[part]]) ** 2, axis=(1, 2)))
        increment_rms[part] = np.sqrt(np.mean((delta[left[part]] - delta[right[part]]) ** 2, axis=(1, 2)))
    _, singular, vh = np.linalg.svd(delta.reshape(-1, f.shape[-1]), full_matrices=False)
    rank = int(np.sum(singular > singular[0] * 1e-12)) if singular.size and singular[0] else 0
    basis = vh[:rank].T
    f0_difference = f[left, anchor] - f[right, anchor]
    projection_norm = np.linalg.norm(f0_difference @ basis, axis=1)
    return dict(pair_left=left, pair_right=right, complete_response_whitened_real_RMS=response_rms,
                target_increment_whitened_real_RMS=increment_rms,
                background_difference_projection_into_full_increment_space_norm=projection_norm,
                full_increment_space_singular_values=singular,
                full_increment_space_empirical_rank=rank)


def offline_distances(prior, response, y, reference, static, selected, consensus_q):
    """All group contrasts include the measured background residual term."""
    ll = reference["log_likelihood"]
    d = ll - ll[:, [prior.central_target]]
    identity = residual_contrast(response, y, prior.central_target)
    left, right = static["pair_left"], static["pair_right"]
    differences = d[left] - d[right]
    tv = 0.5 * np.sum(np.abs(reference["conditionals"][left] - reference["conditionals"][right]), axis=1)
    means = reference["conditionals"] @ prior.pose_nodes
    decision = np.sum((means[left] - means[right]) ** 2, axis=1)
    f = real_white(response)
    delta = f - f[:, [prior.central_target]]
    # Exact two-background identity: <rbar_X,deltaDelta>-<deltaf0,barDelta>.
    identity_pair_error = 0.0
    for start in range(0, len(left), 128):
        part = slice(start, start + 128)
        a, b = left[part], right[part]
        rbar_x = real_white(y) - 0.5 * (f[a] + f[b])
        first = np.sum(rbar_x * (delta[a] - delta[b]), axis=-1)
        second = np.einsum("pm,pxm->px", f[a, prior.central_target] - f[b, prior.central_target],
                           0.5 * (delta[a] + delta[b]))
        identity_pair_error = max(identity_pair_error, float(np.max(np.abs(first-second-differences[part]))))
    selected_mask = np.zeros(len(f), dtype=bool)
    selected_mask[selected] = True
    ss = selected_mask[left] & selected_mask[right]
    radius = 0.5 * np.sum(np.abs(reference["conditionals"] - consensus_q), axis=1)
    return dict(contrast_residual_aware_Linf=np.max(np.abs(differences), axis=1),
                contrast_residual_aware_oscillation=np.ptp(differences, axis=1),
                conditional_target_TV=tv, conditional_mean_squared_distance=decision,
                conditional_TV_to_sample_consensus=radius,
                sampled_pair_mask=ss,
                all_conditional_TV_diameter=float(tv.max()),
                sampled_conditional_TV_diameter=float(tv[ss].max()) if np.any(ss) else 0.0,
                all_conditional_TV_radius_to_sample_consensus=float(radius.max()),
                unqueried_background_posterior_mass=float(reference["background_weights"][~selected_mask].sum()),
                residual_projection_identity_max_abs_error=float(np.max(np.abs(identity-d))),
                two_background_contrast_identity_max_abs_error=identity_pair_error,
                all_conditional_mean_regret_upper=float(np.max(np.sum(
                    (means - consensus_q @ prior.pose_nodes) ** 2, axis=1))))


def algebra_checks():
    """Deterministic wiring checks plus explicit probability counterexamples.

    These tests are not Maxwell evidence.  The accompanying proof note states
    the general Dobrushin argument, which a finite numerical check cannot prove.
    """
    rng = np.random.default_rng(20261029)
    response = (rng.normal(size=(4, 7, 2, 16)) + 1j*rng.normal(size=(4, 7, 2, 16))) * .01
    y = (rng.normal(size=(2, 16)) + 1j*rng.normal(size=(2, 16))) * .01
    ll = potentials(response, y)
    identity_error = float(np.max(np.abs(residual_contrast(response, y, 3) - (ll-ll[:, [3]]))))
    conditional = np.array([[.8,.1,.1],[.2,.7,.1],[.2,.2,.6]])
    beta, beta_tilde = np.array([.2,.3,.5]), np.array([.7,.2,.1])
    diameter, _ = pair_tv(conditional)
    input_tv = float(.5*np.abs(beta-beta_tilde).sum())
    output_tv = float(.5*np.abs(beta@conditional-beta_tilde@conditional).sum())
    epsilon, eta = .01, 1e-6
    table = np.array([[1.,epsilon],[epsilon,1.]])
    scaled = table * np.array([1.,eta])[:,None]
    before, after = table.sum(axis=0)/table.sum(), scaled.sum(axis=0)/scaled.sum()
    contrast_error = float(np.max(np.abs(np.log(table)-np.log(table[:,[0]])
                                        -(np.log(scaled)-np.log(scaled[:,[0]])))))
    tiny = 1e-5
    completion_a = np.tile([tiny,tiny*tiny], (3,1))
    completion_b = completion_a.copy(); completion_b[1:,1] = 1.
    known = np.array([[True,True],[True,False],[True,False]])
    qa, qb = completion_a.sum(axis=0)/completion_a.sum(), completion_b.sum(axis=0)/completion_b.sum()
    if (identity_error > 1e-10 or output_tv > diameter*input_tv+1e-14 or
            contrast_error > 1e-12 or not np.array_equal(completion_a[known], completion_b[known])):
        raise AssertionError("conditional-consensus algebra wiring failed")
    return dict(status="PASS_ALGEBRA_ONLY_NOT_MAXWELL", residual_identity_max_abs_error=identity_error,
                dobrushin=dict(conditional_diameter=diameter, background_weight_TV=input_tv,
                              mixture_TV=output_tv, contraction_upper=diameter*input_tv),
                scale_counterexample=dict(likelihood=table, row_scaled_likelihood=scaled,
                    conditional_contrast_change_max=contrast_error, target_before=before,
                    target_after=after, target_TV=float(.5*np.abs(before-after).sum())),
                unqueried_completion_counterexample=dict(known_mask=known,
                    completion_A=completion_a, completion_B=completion_b,
                    paid_transcript_identical=True, posterior_A=qa, posterior_B=qb,
                    target_TV=float(.5*np.abs(qa-qb).sum()),
                    scope="abstract finite Gaussian-compatible [0,1] likelihood table, not a Maxwell world"))
