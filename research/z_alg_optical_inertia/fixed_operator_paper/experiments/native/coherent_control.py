"""Equal-energy 11-beam native TM control for the fixed old-noise contract."""
import json
import platform

import numpy as np
from scipy.linalg import eigh, eigvalsh
from flint import acb, acb_mat, arb, ctx
from threadpoolctl import threadpool_limits

from energy import clocks, elapsed, endpoints, arb_frame_signs
from source_index import (
    HERE, STUDY, CONTROL_ID, load, ProofGrid, adj, ball, block, diag,
    diagonal_times, gersh, hcat, positive_dimension,
)
from verify import frob, hi, midpoint


def hermitian(a):
    return (a+a.conj().T)/2


def save(path, result):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def interval_from_endpoints(bounds):
    lower, upper = map(arb, bounds)
    return lower.union(upper)


def trace_real(a):
    return sum((a[i, i].real for i in range(a.nrows())), arb(0))


def strict_norm_upper(matrix):
    # An interval square formed as x*x can straddle zero for residual balls.
    return sum((abs(z).upper()**2 for z in matrix.entries()), arb(0)).sqrt().upper()


def construct_sector(g, b, ell, upper, weight):
    endpoints_exact = [(arb(float(l)), arb(float(u))) for l, u in zip(ell, upper, strict=True)]
    mid = [(l+u)/2 for l, u in endpoints_exact]
    products = [l*u for l, u in endpoints_exact]
    mg = diagonal_times(mid, g)
    pg = diagonal_times(products, g)
    pb = diagonal_times(products, b)
    n = len(ell)
    physical = acb(weight)*(diag([1]*n)-mg-adj(mg)+adj(g)*pg)
    coupling = acb(weight)*(diagonal_times(mid, b)-adj(g)*pb)
    constant = acb(weight)*adj(b)*pb
    return physical, coupling, constant


def strict_signs(matrix):
    p, n = positive_dimension(matrix)
    if p+n != matrix.nrows():
        pp, nn = arb_frame_signs(matrix)
        p, n = max(p, pp), max(n, nn)
    return [p, n] if p+n == matrix.nrows() else None


def spectrum_record(a):
    """Every signed gap comes from a verified shifted inertia test."""
    n = a.nrows()
    signs = strict_signs(a)
    ev = eigvalsh(hermitian(midpoint(a)), check_finite=False)
    record = dict(inertia=signs, numeric_minimum=float(ev[0]), numeric_maximum=float(ev[-1]))
    if signs is None:
        return record
    pos, neg = signs
    signed_gap = float(np.min(abs(ev)))
    for factor in (.99, .9, .75, .5, .25, .125):
        shift = signed_gap*factor
        p = positive_dimension(a-acb(shift)*diag([1]*n))[0] if pos else 0
        m = positive_dimension(a+acb(shift)*diag([1]*n))[1] if neg else 0
        if p == pos and m == neg:
            record['strict_spectral_gap_at_least'] = shift
            if pos == n:
                record['minimum_eigenvalue_lower'] = shift
            if neg == n:
                record['maximum_eigenvalue_upper'] = -shift
            break
    return record


def matrices(g, s, b, ell, upper, weight, count_negative, expected_nu):
    """Use one native A solve, then verified receiver-sized Woodbury."""
    t = clocks()
    A, N, C = construct_sector(g, b, ell, upper, weight)
    assembly = elapsed(t)
    t = clocks()
    right = hcat(N, adj(s))
    inverse_times_rhs = A.solve(right)
    ar = b.ncols()
    state = A.nrows()
    ports = s.nrows()
    left = block(inverse_times_rhs, range(state), range(ar))
    right_part = block(inverse_times_rhs, range(state), range(ar, ar+ports))
    source = C-adj(N)*left
    a_solve = elapsed(t)
    t = clocks()
    middle = diag([1]*ports)+s*right_part
    signs = strict_signs(middle)
    if signs is None:
        raise ArithmeticError('Receiver Woodbury block inertia unresolved.')
    nu = count_negative+signs[0]-ports
    if nu != expected_nu:
        raise ArithmeticError('K inertia disagrees with the frozen box receipt.')
    joined_response = inverse_times_rhs-right_part*middle.solve(s*inverse_times_rhs)
    zeros = acb_mat(ar, ports)
    Q = acb_mat(ar+ports, ar+ports,
                 hcat(C, zeros).entries()+hcat(adj(zeros), diag([1]*ports)).entries())
    Q = Q-adj(right)*joined_response
    return source, Q, dict(full_A_assembly=assembly, one_full_A_solve_and_source_schur=a_solve,
                           receiver_woodbury=elapsed(t), woodbury_inertia=signs,
                           certified_K_inertia=[state-nu, nu])


def proposal(source, power, n_beams, radius):
    s, w = (hermitian(midpoint(A)) for A in (source, power))
    values, vectors = eigh(s, w, check_finite=False)
    leading = values[-n_beams:][::-1]
    if np.any(leading <= 0):
        raise ArithmeticError('Floating proposal lacks the required positive generalized directions.')
    spectral = vectors[:, -n_beams:][:, ::-1]*(radius/np.sqrt(leading))[None, :]
    order = sorted(range(power.ncols()), key=lambda j: int(f'{j:05b}'[::-1], 2))[:n_beams]
    planar = np.eye(power.ncols(), dtype=complex)[:, order]
    return dict(coherent=spectral, plane_wave=planar), order, leading


def sharpen_energy(source, power, beams, radius, initial):
    estimates = sorted((float(v.real) for v in power.solve(source).eig(algorithm='approx')),
                       reverse=True)[:beams]
    brackets = []
    for i, estimate in enumerate(estimates):
        for relative in (1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-6):
            lower, upper = estimate*(1-relative), estimate*(1+relative)
            if lower <= 0:
                break
            p, _ = arb_frame_signs(source-acb(lower)*power)
            _, q = arb_frame_signs(source-acb(upper)*power)
            if p >= i+1 and q >= source.nrows()-i:
                brackets.append([lower, upper])
                break
    if len(brackets) != beams:
        return initial, [], 'previously verified spectral roots'
    bnd = arb(radius)**2*sum((1/arb(upper).union(arb(lower)) for lower, upper in brackets), arb(0))
    tightened = endpoints(bnd)
    if not interval_from_endpoints(initial).overlaps(bnd):
        raise ArithmeticError('Independent energy brackets disagree.')
    return tightened, brackets, 'independently reverified 128-bit shifted pencil signs'


def normalization_scale(base, power, target):
    energy = trace_real(adj(base)*power*base)
    if not energy > 0:
        raise ArithmeticError('Incident fluence of proposed directions unresolved.')
    return (target/energy).sqrt()


def physical_fluence(full_incident, beam0, V):
    return frob(full_incident*V)**2/beam0


def verify_arm(name, direction, source, Q, S, B, G, material, dictionary_current,
               power, full_incident, beam0, energy_interval, absolute_noise,
               perturbation, radius, arrays):
    tick = clocks()
    n, s = direction.shape
    ports = S.nrows()
    base = ball(direction)
    scale = normalization_scale(base, power, energy_interval)
    V = acb(scale)*base
    beam_fluences = acb(scale**2)*(adj(base)*power*base)
    energy = trace_real(beam_fluences)
    direct_energy = scale**2*physical_fluence(full_incident, beam0, base)
    if not direct_energy.overlaps(energy):
        raise ArithmeticError(f'{name} direct full-domain fluence disagrees with the Gram calculation.')
    Vnom = ball(midpoint(V))
    omega = acb(arb(radius)**2)*diag([1]*s)
    source_test = acb(scale**2)*(adj(base)*source*base)-omega
    source_record = spectrum_record(source_test)
    if source_record['inertia'] is None:
        raise ArithmeticError(f'{name} source sign check unresolved.')

    # One verified dictionary current is reused by both coherent and planar
    # mixtures. Nominal points are stored alongside their rigorous envelopes.
    J = acb(scale)*(dictionary_current*base)
    clean = S*J
    ynom = midpoint(clean)
    measured = ynom+perturbation
    measured_ball = ball(measured)
    error = measured_ball-clean
    column_bounds = [frob(block(error, range(ports), [i])) for i in range(s)]
    if not all(r <= absolute_noise for r in column_bounds) or not frob(error) <= arb(radius):
        raise ArithmeticError(f'{name} violates the fixed measurement error contract.')
    Jnom = ball(midpoint(J))
    incident = B*Vnom
    constitutive_defect = Jnom-diagonal_times(material, incident+G*Jnom)

    aa = block(Q, range(n), range(n))
    ay = block(Q, range(n), range(n,n+ports))
    ya = block(Q, range(n,n+ports), range(n))
    yy = block(Q, range(n,n+ports), range(n,n+ports))
    cross = adj(base)*ay*measured_ball
    old_form = (acb(scale**2)*(adj(base)*aa*base)+acb(scale)*(cross+adj(cross))
                +adj(measured_ball)*yy*measured_ball-omega)
    po_record = spectrum_record(old_form)
    if po_record['inertia'] != [0, s]:
        raise ArithmeticError(f'{name} nonempty data did not yield certified zero-index saturation.')
    mix = acb(scale)*(ya*base)+yy*measured_ball
    Ry = yy-mix*old_form.solve(adj(mix))
    receiver_record = spectrum_record(Ry)
    source_deficit = s-source_record['inertia'][0]
    positive_theorem = ports-source_deficit
    receiver_record['verified_positive_by_index_identity'] = positive_theorem
    if receiver_record['inertia'] is not None and receiver_record['inertia'][0] != positive_theorem:
        raise ArithmeticError(f'{name} direct conditioned receiver count disagrees with index formula.')

    arrays[f'{name}_directions'] = direction
    arrays[f'{name}_V_midpoint'] = midpoint(V)
    arrays[f'{name}_J_midpoint'] = midpoint(Jnom)
    arrays[f'{name}_Y_clean_midpoint'] = ynom
    arrays[f'{name}_Y_measured'] = measured
    arrays[f'{name}_source_threshold_matrix_midpoint'] = midpoint(source_test)
    arrays[f'{name}_old_threshold_matrix_midpoint'] = midpoint(old_form)
    arrays[f'{name}_Ry_midpoint'] = midpoint(Ry)
    return dict(
        arm=name, number_of_beams=s, target_exact_equal_by_normalization=True,
        fluence_interval=endpoints(energy), nominal_fluence_interval=endpoints(trace_real(adj(Vnom)*power*Vnom)),
        direct_full_grid_fluence_interval=endpoints(direct_energy),
        per_beam_fluence_intervals=[endpoints(beam_fluences[i,i].real) for i in range(s)],
        source_deficit_at_this_acquisition=source_deficit,
        source_test=source_record, old_data_threshold=po_record, receiving_block=receiver_record,
        prediction_status='BOUNDED_RECEIVER_ELLIPSOID' if positive_theorem == ports
                          else 'POSITIVE_RECEIVER_BLOCK_OBSTRUCTED',
        admissible_same_material=True, exact_forward_solution_from_verified_dictionary_solve=True,
        physical_nominal_current_residual_frobenius_upper=hi(strict_norm_upper(constitutive_defect)),
        physical_nominal_current_relative_residual_upper=hi(strict_norm_upper(constitutive_defect)/frob(diagonal_times(material, incident))),
        worst_measured_error_column_upper=max(hi(r) for r in column_bounds),
        measured_error_frobenius_upper=hi(frob(error)),
        stored_clean_nominal_error_upper=hi(strict_norm_upper(ball(ynom)-clean)),
        timing=elapsed(tick),
    )


def occupied_wide_control(energy, disk, cells, profile, g, s, b, dictionary_current):
    started = clocks()
    r0 = next(r for r in energy['rows'] if r['box'] == 'solid_half_2.5')
    width, weight = r0['half_width_or_cap'], r0['weight']
    q, dictionary_positive = r0['negative_material_modes'], r0['positive_dictionary_modes']
    ell, upper = np.maximum(0., profile-width), profile+width
    if (q, dictionary_positive) != (13, 6) or not np.all(ell <= profile) or not np.all(profile <= upper):
        raise AssertionError('Wide occupied control contract changed.')
    source, Q, costs = matrices(g, s, b, ell, upper, weight, q, expected_nu=13)
    m, ports = b.ncols(), s.nrows()
    nu = costs['certified_K_inertia'][1]
    radius = disk.radius
    noise_seed = 202611261013
    rng = np.random.default_rng(noise_seed)
    noise = rng.normal(size=(ports, m))+1j*rng.normal(size=(ports, m))
    noise *= (.9*radius/np.sqrt(m)/np.linalg.norm(noise, axis=0))[None,:]
    order = sorted(range(m), key=lambda j: int(f'{j:05b}'[::-1],2))
    clean = s*block(dictionary_current, range(len(cells)), order)
    clean_center = midpoint(clean)
    y = clean_center+noise
    actual_error = ball(y)-clean
    radius_column = arb(radius)/arb(m).sqrt()
    errors = [strict_norm_upper(block(actual_error, range(ports), [i])) for i in range(m)]
    if not all(r <= radius_column for r in errors):
        raise ArithmeticError('Wide occupied physical measurements exceed column noise caps.')
    output = dict(
        precision_bits=ctx.prec, status='COMPLETE', geometry=[96,96], supported_cells=len(cells),
        box='solid_half_2.5', weight=weight, material='same stored c_solid used by the narrow control',
        common_physical_solve='coherent_control.json dictionary_forward, 32 RHS; no additional physical solves',
        source_only_dictionary_positive=dictionary_positive, material_negative=q, state_negative=nu,
        radius=radius, noise_seed=noise_seed, maximum_verified_noise_column=hi(max(errors)),
        noise_contract='r^2 I_s on both nested prefixes, with 32 absolute column caps r/sqrt(32)',
        illumination='unit-amplitude bit-reversal prefixes from the fixed dictionary', rows=[], costs=costs,
    )
    for count in (16, 32):
        tick = clocks()
        selected = order[:count]
        V = np.eye(m)[:,selected]
        old = ball(np.vstack((V, y[:,:count])))
        omega = acb(arb(radius)**2)*diag([1]*count)
        src = spectrum_record(block(source, selected, selected)-omega)
        Po = adj(old)*Q*old-omega
        old_proof = spectrum_record(Po)
        record = dict(s=count, source_threshold=src, old_threshold=old_proof,
                      verified_error_frobenius_upper=hi(strict_norm_upper(
                          block(actual_error, range(ports), range(count)))),
                      dictionary_shortfall_lower=q-dictionary_positive)
        if src['inertia'] is not None and old_proof['inertia'] is not None:
            h = q-src['inertia'][0]
            k = old_proof['inertia'][0]
            if h < q-dictionary_positive or k > nu:
                raise ArithmeticError('Wide measurement contradicts its proved nonempty contract.')
            mix = block(Q, range(m,m+ports), range(m+ports))*old
            Ry = block(Q, range(m,m+ports), range(m,m+ports))-mix*Po.solve(adj(mix))
            ry = spectrum_record(Ry)
            p = ports-h+nu-k
            if ry['inertia'] is not None and ry['inertia'][0] != p:
                raise ArithmeticError('Wide direct receiver count differs from the index identity.')
            record.update(k=k, h=h, saturated=k == nu, receive=ry,
                          verified_positive_by_index_identity=p,
                          bounded_receiver_ellipsoid=(k == nu and p == ports),
                          status='SATURATED' if k == nu else 'UNSATURATED')
        else:
            record['status'] = 'UNRESOLVED_SPECTRAL_SIGNS'
        record['timing'] = elapsed(tick)
        output['rows'].append(record)
    tick = clocks()
    np.savez_compressed(HERE/'occupied_wide_control_arrays.npz', cells=cells, material=profile,
                        Y_clean_midpoint=clean_center, Y_measured=y, source_coefficients=np.eye(m)[:,order],
                        original_dictionary_indices=order)
    costs['arrays_persistence'] = elapsed(tick)
    costs['total_incremental'] = elapsed(started)
    save(HERE/'occupied_wide_control.json', output)
    print('occupied wide', [(r['s'], r.get('k'), r.get('h'), r.get('verified_positive_by_index_identity'))
                            for r in output['rows']], costs['total_incremental'], flush=True)


def main():
    started = clocks()
    energy = json.loads((HERE/'energy.json').read_text())
    row = next(r for r in energy['rows'] if r['box'] == 'solid_half_1.5')
    disk = load(STUDY/'ip_guided_v4/boxes/results', [CONTROL_ID])[CONTROL_ID]
    with np.load(STUDY/'ip_final_v5/strong_control/strong_ring_arrays.npz', allow_pickle=False) as arrays:
        cells, midpoint_material = arrays['ids'], arrays['c_solid']
    weight, width, beams = row['weight'], row['half_width_or_cap'], row['negative_material_modes']
    if (len(cells), weight, beams) != (1153, 2**-13, 11) or disk.geometry.side != 96:
        raise AssertionError('Native narrow-box control changed.')
    radius = disk.radius
    inf_bounds = row['energy']['disk']['interval']
    per_beam_radius = arb(radius)/arb(beams).sqrt()
    ell = np.maximum(0., midpoint_material-width)
    upper = midpoint_material+width
    if not np.all(ell <= midpoint_material) or not np.all(midpoint_material <= upper):
        raise AssertionError('Generator is outside the material box.')
    _, (g, s, b) = ProofGrid(disk.geometry).restricted(cells)
    bd = ball(disk.geometry.B)
    power = adj(bd)*bd
    beam0 = power[0,0].real
    power = power/acb(beam0)
    if strict_signs(power) != [b.ncols(), 0]:
        raise ArithmeticError('Incident full-domain fluence Gram not proved positive.')
    common_setup = elapsed(started)

    source, Q, matrix_costs = matrices(g, s, b, ell, upper, weight, beams, expected_nu=0)
    tick = clocks()
    sharpened, roots, energy_method = sharpen_energy(source, power, beams, radius, inf_bounds)
    target = arb(5)*interval_from_endpoints(sharpened)/4
    proposals, order, generalized = proposal(source, power, beams, radius)
    seed = 202611260959
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=(s.nrows(), beams))+1j*rng.normal(size=(s.nrows(), beams))
    noise *= (0.9*radius/np.sqrt(beams)/np.linalg.norm(noise, axis=0))[None,:]
    design_cost = elapsed(tick)

    output = dict(
        precision_bits=ctx.prec, arb_threads=ctx.threads, blas_threads=2, platform=platform.machine(),
        status='DESIGNS_FROZEN', box='solid_half_1.5', supported_cells=len(cells), grid=[96,96],
        material='stored c_solid midpoint profile, exact exterior vacuum',
        weight=weight, full_dictionary_size=b.ncols(), n_receivers=s.nrows(),
        error_ball_radius=radius, exact_per_beam_noise_cap='r / sqrt(11)',
        numeric_per_beam_noise_cap_enclosure=endpoints(per_beam_radius),
        omega='r^2 I_11 for both arms, justified by 11 column caps',
        realized_noise='shared 0.9 nominal cap per column, verified against exact full-wave responses',
        common_exact_total_fluence_definition='(5/4) times exact E_inf for the disk noise radius',
        disk_e_inf_interval=sharpened, previous_disk_e_inf_interval=inf_bounds,
        independently_verified_energy_roots=roots, energy_reverification=energy_method,
        common_target_fluence_interval=endpoints(target),
        normalizing_metric='trace(V* W_in V), W_in=B_full* B_full / column-0 incident fluence',
        plane_wave_zero_based_dictionary_indices=order, numeric_generalized_eigenbeam_proposals=generalized.tolist(),
        numerical_operator_contract='EXACT_STORED_FLOAT of padded spectrum and source/receiver arrays',
        source_rank_provenance='strict n_minus(A)=11 in energy.json for this same box, sector and stored operator',
        cost=dict(setup=common_setup, **matrix_costs, freeze_designs=design_cost), arms=[],
    )
    files = dict(cells=cells, material_in_roi=midpoint_material, verified_operator_matrix_G_midpoint=midpoint(g),
                 incident_dictionary_in_roi=disk.geometry.B[cells], receiver_in_roi=disk.geometry.S[:, cells],
                 power_gram_midpoint=midpoint(power), common_noise_proposal=noise,
                 source_pencil_midpoint=midpoint(source), full_port_form_midpoint=midpoint(Q))
    for name, directions in proposals.items():
        files[f'{name}_directions'] = directions
    np.savez_compressed(HERE/'coherent_control_design.npz',
                        **{f'{name}_directions':v for name,v in proposals.items()},
                        plane_wave_indices=order, radius=radius, target_fluence_bounds=endpoints(target))
    save(HERE/'coherent_control.json', output)
    print('frozen source designs', output['common_target_fluence_interval'], matrix_costs, flush=True)

    tick = clocks()
    chi = [arb(float(v)) for v in midpoint_material]
    dictionary_current = (diag([1]*len(cells))-diagonal_times(chi, g)).solve(diagonal_times(chi, b))
    dict_rhs = diagonal_times(chi, b)
    dict_check = dictionary_current-diagonal_times(chi, b+g*dictionary_current)
    output['dictionary_forward'] = dict(
        one_full_solve=True, simultaneous_right_hand_sides=dictionary_current.ncols(),
        invertible_by_verified_solve=True, ambient_exterior_equations_exact_zero=True,
        enclosure_of_original_equation_residual_upper=hi(strict_norm_upper(dict_check)),
        enclosure_of_original_equation_relative_residual_upper=hi(strict_norm_upper(dict_check)/frob(dict_rhs)),
        timing=elapsed(tick),
    )
    print('verified physical midpoint dictionary', output['dictionary_forward'], flush=True)

    for name, direction in proposals.items():
        outcome = verify_arm(name, direction, source, Q, s, b, g, chi, dictionary_current,
                             power, bd, beam0, target, per_beam_radius, noise, radius, files)
        output['arms'].append(outcome)
        save(HERE/'coherent_control.json', output)
        print(name, outcome['source_test'], outcome['receiving_block'], flush=True)
    tick = clocks()
    files['dictionary_full_wave_currents_midpoint'] = midpoint(dictionary_current)
    np.savez_compressed(HERE/'coherent_control_arrays.npz', **files)
    output['status'] = 'CERTIFIED'
    output['cost']['arrays_persistence'] = elapsed(tick)
    output['cost']['total'] = elapsed(started)
    output['noise_proposal_seed'] = seed
    save(HERE/'coherent_control.json', output)
    print('certified', output['cost']['total'], flush=True)
    occupied_wide_control(energy, disk, cells, midpoint_material, g, s, b, dictionary_current)


if __name__ == '__main__':
    ctx.prec, ctx.threads = 128, 3
    with threadpool_limits(limits=2):
        main()
