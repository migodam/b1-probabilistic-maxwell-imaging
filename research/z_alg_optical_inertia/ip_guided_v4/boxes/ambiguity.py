"""Independently verify two archived feasible Maxwell materials in disjoint boxes."""
import argparse
import importlib.util
from pathlib import Path
import platform

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

from run import HERE, STUDY, clocks, read, timing, write
from verify import ball, diagonal_times, frob, frob2, hi, lo, midpoint

CASE = 'GENERAL6_00'
WORKSPACE = STUDY.parents[1]
ORIGINAL = STUDY.parent / 'delegated/z_alg_r2_publication/repository/research/research_z_alg_r2'
OBSERVATION = Path('experiments/general6_public') / (CASE + '.npz')
GEOMETRY = Path('experiments/GENERAL6_CONTRACT.json')
WITNESSES = Path('results') / ('general6_' + CASE + '_general_quotient20.json')
SOURCE = Path('src/physics.py')


def archived_geometry(root):
    contract = read(root / GEOMETRY)
    spec = importlib.util.spec_from_file_location('general6_r2_physics', root / SOURCE)
    physics = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(physics)
    model = physics.Model(**contract['geometry'])
    with np.load(root / OBSERVATION, allow_pickle=False) as archive:
        Y = archive['y']
    data = read(root / WITNESSES)
    if data['case'] != CASE or data['method'] != 'general_quotient20':
        raise ArithmeticError('The archived root table does not belong to the selected observation.')
    rows = sorted((row for row in data['independently_validated'] if row['status'] == 'QUALIFIED'
                   and row['validation']['qualified']), key=lambda v: v['source_path'])
    if len(rows) != 2 or [r['source_path'] for r in rows] != [16, 17]:
        raise ArithmeticError('The original Q20 paths 16 and 17 are required.')
    provenance = {key: str((root / path).relative_to(WORKSPACE)) if
                  (root / path).is_relative_to(WORKSPACE) else str(root / path)
                  for key, path in (('observation', OBSERVATION), ('geometry_contract', GEOMETRY),
                                    ('archived_candidates', WITNESSES), ('operator_source', SOURCE))}
    provenance.update(observation_key='y', source_path_ids=[r['source_path'] for r in rows],
                      geometry_schema=contract['schema'], geometry=contract['geometry'],
                      originally_noiseless=contract['observations']['noise'] == 'none',
                      coefficient_platform=platform.platform(),
                      exact_coefficients='Binary64 G, S and E0 produced by the unmodified archived pulse-TM Model')
    return dict(G=model.G, S=model.S, B=model.E0, Y=Y,
                center_materials=np.asarray([r['m'] for r in rows], dtype=float),
                source_path_ids=np.asarray([r['source_path'] for r in rows], dtype=int),
                material_box=np.asarray(contract['material_box'], dtype=float)), provenance


def outward_boxes(materials, allowed):
    half_width = arb(1)/20
    centers = [[arb(float(x)) for x in vector] for vector in materials]
    low = [[x-half_width for x in vector] for vector in centers]
    high = [[x+half_width for x in vector] for vector in centers]
    if not all(low_i >= arb(float(allowed[0])) and high_i <= arb(float(allowed[1]))
               for lower, upper in zip(low, high, strict=True)
               for low_i, high_i in zip(lower, upper, strict=True)):
        raise ArithmeticError('An exact ±0.05 box extends outside the original admissible material range.')
    displacement = [abs(a-b) for a, b in zip(*centers, strict=True)]
    gaps = [max(arb(0), (d-2*half_width).lower()) for d in displacement]
    largest = max(gaps)
    if not largest > 0:
        raise ArithmeticError('The exact boxes are not strictly disjoint.')
    l2 = sum((v*v for v in displacement), arb(0)).sqrt()
    gap_l2 = sum((v*v for v in gaps), arb(0)).sqrt().lower()
    rms = l2/arb(len(displacement)).sqrt()
    rms_gap = gap_l2/arb(len(displacement)).sqrt()
    separated_axis = int(np.argmax([float(v) for v in gaps]))
    pair = dict(exact_box_half_width='1/20', box_half_width=0.05,
                box_centers_admissible=True, both_complete_boxes_inside_original_material_range=True,
                boxes_disjoint=True, max_separating_coordinate_zero_based=separated_axis,
                max_coordinate_box_gap_lower=lo(largest),
                center_ell2_distance_interval=[lo(l2), hi(l2)],
                center_rms_distance_interval=[lo(rms), hi(rms)],
                center_ellinf_distance_interval=[lo(max(displacement)), hi(max(displacement))],
                minimum_boxes_ell2_distance_lower=lo(gap_l2),
                minimum_boxes_rms_distance_lower=lo(rms_gap),
                minimax_single_estimator_rms_lower=lo(rms/2))
    arrays = dict(box_lower_outward=np.asarray([[lo(x) for x in v] for v in low]),
                  box_upper_outward=np.asarray([[hi(x) for x in v] for v in high]))
    return centers, pair, arrays


def physical_check(inputs):
    G, S, B, Y = (ball(inputs[k]) for k in ('G', 'S', 'B', 'Y'))
    material = inputs['center_materials']
    centers, pair, arrays = outward_boxes(material, inputs['material_box'])
    Y_norm = frob(Y)
    radius = Y_norm / 2**30
    if not radius > 0:
        raise ArithmeticError('A nonzero observation is required for the relative noise ball.')
    I = ball(np.eye(G.nrows()))
    proof, fields, currents, ports, exact_ports = [], [], [], [], []
    for index, (m, archived_path) in enumerate(zip(centers, inputs['source_path_ids'], strict=True)):
        started = clocks()
        A = I-diagonal_times(m, G)
        J = A.solve(diagonal_times(m, B))
        E = B + G*J
        P = S*J
        residual = frob(P-Y)
        distance = residual / radius
        if not residual.upper() <= radius.lower():
            raise ArithmeticError(f'Archived path {archived_path} failed the chosen original-data ball.')
        state_error = frob2(A*J-diagonal_times(m, B)).upper().sqrt()
        proof.append(dict(box=f'box_{index+1}', source_path=int(archived_path), center=material[index].tolist(),
                          full_state_matrix_verified_invertible=True, center_feasible=True,
                          interval_rhs=B.ncols(), full_state_dimension=G.nrows(),
                          residual_absolute_lower=max(0., lo(residual)),
                          residual_absolute_upper=hi(residual), residual_over_radius_lower=max(0., lo(distance)),
                          residual_over_radius_upper=hi(distance), state_defect_frobenius_upper=hi(state_error),
                          **timing(started)))
        currents.append(midpoint(J))
        fields.append(midpoint(E))
        ports.append(midpoint(P))
        exact_ports.append(P)
    ports_difference = frob(exact_ports[0]-exact_ports[1])
    pair.update(center_ports_difference_absolute_interval=[max(0., lo(ports_difference)), hi(ports_difference)],
                both_boxes_have_feasible_materials=True)
    arrays.update(current_matrices=np.stack(currents), total_field_matrices=np.stack(fields),
                  predicted_ports=np.stack(ports), radius_enclosure=np.array([lo(radius), hi(radius)]))
    return dict(observed_Y_norm_interval=[lo(Y_norm), hi(Y_norm)], radius_absolute_interval=[lo(radius), hi(radius)],
                radius_relative_to_Y=2.**-30, radius_definition='2^-30 times exact Frobenius norm of archived Y',
                witnesses=proof, material_separation=pair, full_forward_interval_solves=len(centers),
                total_forward_interval_rhs=sum(v['interval_rhs'] for v in proof)), arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=ORIGINAL)
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--replay-npz', type=Path, help='replay the embedded exact represented inputs')
    args = parser.parse_args()
    ctx.prec, ctx.threads = 128, 1
    start = clocks()
    if args.replay_npz is None:
        inputs, provenance = archived_geometry(args.archive.resolve())
        provenance_mode = 'archived geometry and observations'
    else:
        with np.load(args.replay_npz, allow_pickle=False) as f:
            inputs = {k: f[k] for k in ('G', 'S', 'B', 'Y', 'center_materials', 'source_path_ids', 'material_box')}
        provenance = read(args.replay_npz.with_suffix('.json'))['provenance']
        provenance_mode = 'embedded exact binary64 arrays'
    verdict, arrays = physical_check(inputs)
    result = dict(observation_id=CASE, material_count=int(inputs['G'].shape[0]),
                  receiver_count=int(inputs['S'].shape[0]), source_count=int(inputs['B'].shape[1]),
                  arith_precision_bits=ctx.prec, provenance=provenance, provenance_mode=provenance_mode,
                  status='two_disjoint_boxes_each_contain_a_verified_feasible_center',
                  scope='Existence of two distinct discrete whole-wave material explanations in one shared positive-radius ball.',
                  **verdict, **timing(start))
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out / 'ambiguity.npz', **inputs, **arrays)
    write(args.out / 'ambiguity.json', result)
    print(CASE, 'residual/r', [r['residual_over_radius_upper'] for r in result['witnesses']],
          'pair RMS', result['material_separation']['center_rms_distance_interval'], flush=True)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
