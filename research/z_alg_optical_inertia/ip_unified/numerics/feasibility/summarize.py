"""Collect auditable physical upper endpoints and failure classifications for all requested boxes."""
import csv
import json
from pathlib import Path

import numpy as np
from flint import ctx

from enclose_native import CONTROL_ID, ball, hi, norm
from search_native import CAPS, HERE, STUDY, inputs, stem, write_json


def read(p):
    return json.loads(p.read_text())


def load_paths(names):
    return {p.name: [read(f) for f in sorted(p.glob('optfresh*__u*.json'))]
            for p in names if p.is_dir()}


def main():
    ctx.prec, ctx.threads = 128, 1
    _, cases = inputs()
    selection = load_paths([HERE/'search_runs', HERE/'projected_runs', HERE/'cie_runs'])
    proof = load_paths([HERE/'strict_endpoints', HERE/'strict_projected'])
    control = HERE/'strict_endpoints'/f'{CONTROL_ID}.json'
    rows = []
    for obs in cases.values():
        vacuum_distance_upper = hi(norm(ball(obs.y)))
        for cap in CAPS:
            same_cap_attempts = [dict(directory=directory, **r) for directory, results in selection.items() for r in results
                                 if r['observation'] == obs.id and r['cap'] == cap]
            all_strict_options = [dict(directory=directory, **r) for directory, results in proof.items() for r in results
                                  if r['observation'] == obs.id and r['cap'] <= cap and 'distance_upper_for_fixed_candidate' in r]
            strict = min(all_strict_options, key=lambda x: x['distance_upper_for_fixed_candidate'], default=None)
            numeric = min(same_cap_attempts, key=lambda x: x['numerical_residual'], default=None)
            if strict and strict['distance_upper_for_fixed_candidate'] <= vacuum_distance_upper:
                bound = strict['distance_upper_for_fixed_candidate']
                source = f"{strict['directory']}/{stem(obs, strict['cap'])}.json"
            else:
                bound = vacuum_distance_upper
                source = 'exact vacuum field J=0'
            proof_pending = numeric is not None and numeric['numerical_residual'] < obs.padded_radius
            classification = ('certified_physical_feasible_witness' if bound <= obs.padded_radius else
                              'numerical_interior_pending_enclosure' if proof_pending else
                              'tested_trajectories_no_feasible_witness' if same_cap_attempts else
                              'no_completed_trajectory')
            rows.append(dict(observation=obs.id, cap=cap, radius=obs.padded_radius, classification=classification,
                             absolute_upper_bound_on_best_residual=bound, upper_over_radius=bound/obs.padded_radius,
                             strict_upper_bound_source=source,
                             best_observed_numerical_residual=None if numeric is None else numeric['numerical_residual'],
                             best_observed_numerical_residual_over_radius=None if numeric is None else numeric['residual_over_radius'],
                             best_numeric_source=None if numeric is None else f"{numeric['directory']}/{stem(obs,cap)}.json",
                             completed_searches=len(same_cap_attempts)))
    offline = 'physical_revision/array/forward_audit.json'
    historical_nonvacuity = [row for row in read(STUDY/offline)['rows'] if row['id'] in {o.id for o in cases.values()}]
    failed = HERE/'initial_vacuum_audit_attempts.json'
    bad = read(failed) if failed.exists() else []
    cie_setup = [read(f) for f in sorted((HERE/'cie_runs').glob('*__cie_setup.json'))]
    projected_stop = HERE/'projected_stop_accounting.txt'
    projected_setup = HERE/'projected_runs'/'setup.json'
    audit_log = HERE/'strict_endpoints'/'audit_invocations.jsonl'
    report = dict(geometry='400 MHz, 2 m square, 96² full-domain TM, equally normalized 32tx and 32rx',
                  observations='Historical paired full 2% Frobenius balls; saved Y only for the online search',
                  actual_model_assumption='Exact represented native Fourier operator, real unchanged box',
                  physical_lower_endpoint_from_search=None,
                  cie_shared_proposal_setup=cie_setup,
                  projected_shared_setup=read(projected_setup) if projected_setup.exists() else None,
                  projected_process_total_including_end_interruption=projected_stop.read_text().strip()
                      if projected_stop.exists() else None,
                  costs_includes_completed_searches_and_completed_proofs_only=True,
                  search_costs_by_directory={name: dict(cases=len(records), cpu_seconds=sum(r['cpu_seconds'] for r in records),
                     summed_case_wall_seconds=sum(r['wall_seconds'] for r in records),
                     peak_rss_bytes=max((r['peak_rss_bytes'] for r in records), default=None),
                     native_physical_calls=sum(r['counts']['native_full_wave_calls'] for r in records),
                     native_fft_columns=sum(r['counts']['fft_columns'] for r in records),
                     state_rhs=sum(r['counts'].get('state_rhs',0) for r in records),
                     adjoint_rhs=sum(r['counts'].get('adjoint_rhs',0) for r in records)) for name, records in selection.items()},
                  logged_strict_audit_invocations=[json.loads(x) for x in audit_log.read_text().splitlines()]
                      if audit_log.exists() else [],
                  proof_costs_by_directory={name: dict(cases=len(records), cpu_seconds=sum(r['cpu_seconds'] for r in records),
                      summed_case_wall_seconds=sum(r['wall_seconds'] for r in records),
                      interval_green_columns=sum(r.get('forward', {}).get('interval_green_columns',0) for r in records))
                      for name, records in proof.items()},
                  rejected_initial_vacuum_audits=bad,
                  historical_generating_controls_source=offline,
                  historical_generating_controls=historical_nonvacuity,
                  separate_low_peak_native_control=read(control) if control.exists() else None,
                  rows=rows)
    write_json(HERE/'upper_endpoints.json', report)
    with (HERE/'upper_endpoints.csv').open('w', newline='') as f:
        out = csv.DictWriter(f, fieldnames=list(rows[0]))
        out.writeheader()
        out.writerows(rows)


if __name__ == '__main__':
    main()
