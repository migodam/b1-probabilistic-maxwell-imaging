"""Explicit recovery of a pre-Book parser failure; identical paid job IDs.

The original failed sequence, source and stderr remain untouched. This does
not reset authorization, ledger or paid jobs; only one metadata keyword was
renamed so Python can parse the runner. No scientific setting changed.
"""
import json
import run_p3l_surrogate_sequence as sequence


def main():
    original=sequence.FOLDER
    failed=json.loads((original/'CLOSED.json').read_text())
    if failed['status']!='ERROR':raise RuntimeError('explicit retained parser failure required')
    if any(any(cost[k] for k in ('RHS','states','LU')) for cost in failed['all_planned_job_costs_including_failures'].values()):
        raise RuntimeError('physical cost already exists; this parser-only recovery is not applicable')
    if any((sequence.OUT/'C/jobs'/task['job']).exists() for task in sequence.tasks()):
        raise RuntimeError('paid/partial job exists; no replay allowed')
    sequence.FOLDER=sequence.OUT/'schedules/surrogate_force_hmc_fix1'
    sequence.OWNED_SOURCE=sequence.OWNED_SOURCE+['tools/run_p3l_surrogate_sequence_fix1.py']
    sequence.main()


if __name__=='__main__':main()
