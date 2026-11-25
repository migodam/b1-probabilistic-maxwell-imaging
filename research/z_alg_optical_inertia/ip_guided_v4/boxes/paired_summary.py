"""Extract the strict paired box rows in observation order for the manuscript."""
import argparse
import csv
from pathlib import Path

from run import CONTROL_ID, HERE, LEGACY_ID, NATIVE_IDS, SMALL_ID, read


FIELDS = ('observation_id', 'grid_side', 'hypothesis', 'observation_radius',
          'active_cells', 'unknown_cells', 'fixed_nonzero_cells', 'support_radius_fraction',
          'center_x_fraction', 'center_y_fraction', 'box_half_width',
          'peak_lower', 'peak_upper', 'T_norm_lower', 'T_spectral_radius_lower',
          'negative_directions', 'sources', 'sources_retained', 'rank_tail_over_radius_lower',
          'positive_same_box_over_radius_lower', 'positive_same_box_scope',
          'zero_lower_same_support_over_radius_lower', 'other_same_support_physical_fit_over_radius_upper',
          'strict_residual_center_over_radius_upper', 'locked_multiplier', 'endpoint_cpu_seconds')


def row(case, rank, *, hypothesis, host_cells=0, positive=None, broad=None, physical_fit=None, rho_lower=None):
    params = rank.get('parameters', rank)
    coordinates = params.get('placement_center', params.get('center'))
    if rho_lower is None:
        rho_lower = (rank.get('homogeneous_endpoint_moment') or {}).get('rho_lower')
    dark = rank.get('exact_dark_check') or {}
    positive_scope = ('obstructed for every scalar mu>0 at W=I with sector alone'
                      if dark.get('every_mu_positive_indefinite_as_lower_direction')
                      else ('same sampled box and state weight' if positive else None))
    return dict(grid_side=case['side'], hypothesis=hypothesis, observation_radius=case['radius'],
                active_cells=rank['cells'], unknown_cells=rank['cells']-host_cells, fixed_nonzero_cells=host_cells,
                support_radius_fraction=params['radius_fraction'], center_x_fraction=coordinates[0],
                center_y_fraction=coordinates[1], box_half_width=params['width'],
                peak_lower=rank['peak_lower'], peak_upper=rank['peak_upper'],
                T_norm_lower=rank['center']['T_norm_lower'], T_spectral_radius_lower=rho_lower,
                negative_directions=rank['nu'], sources=rank['source_dimension'],
                sources_retained=rank['source_dimension']-rank['nu'],
                rank_tail_over_radius_lower=rank['distance_over_radius_lower'],
                positive_same_box_over_radius_lower=positive['distance_over_radius_lower'] if positive else None,
                positive_same_box_scope=positive_scope,
                zero_lower_same_support_over_radius_lower=broad['distance_over_radius_lower'] if broad else None,
                other_same_support_physical_fit_over_radius_upper=physical_fit,
                strict_residual_center_over_radius_upper=rank['center']['center_distance_over_radius_upper'],
                locked_multiplier=rank['mu'], endpoint_cpu_seconds=rank['cost']['cpu_seconds'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=HERE / 'results')
    args = ap.parse_args()
    res = args.out
    primary = read(res / 'verified.json')['cases']
    enlarged = read(res / 'followup_broad.json')['cases']
    low = read(res / 'followup_disk.json')['cases'][CONTROL_ID]['strict_checks']['rank_tail_disk']
    supplemental = read(res / 'followup_host.json')['cases']
    disk_radius = read(res / 'disk_radius.json')
    records = []
    for name in (LEGACY_ID, SMALL_ID, *NATIVE_IDS):
        group = primary[name]['tests']
        broad = enlarged[name]['strict_checks']['broad']
        records.append(dict(observation_id=name, **row(primary[name], group['rank_tail'],
                                                       hypothesis='vacuum fixed outside support',
                                                       positive=group['same_box_positive'], broad=broad)))
    records.append(dict(observation_id=CONTROL_ID,
                        **row(primary[CONTROL_ID], low, hypothesis='disk control, shifted high-peak box',
                              physical_fit=low['same_support_low_physical_witness']['residual_over_radius_upper'],
                              rho_lower=disk_radius.get('spectral_check', {}).get(
                                  'ideal_exclusion_width_operator_rho_lower'))))
    for name in NATIVE_IDS:
        for suffix in ('identity', 'inverse_center_illumination'):
            rank = supplemental[name]['strict_checks']['rank_tail_' + suffix]
            records.append(dict(observation_id=name, **row(primary[name], rank,
                                                           hypothesis='fixed nonzero host, ' + suffix,
                                                           host_cells=rank['parameters']['host_cells'])))
    with (res / 'paired_summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, FIELDS)
        writer.writeheader()
        writer.writerows(records)


if __name__ == '__main__':
    main()
