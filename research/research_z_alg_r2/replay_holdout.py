"""Replay the frozen four-arm benchmark into a new directory, without reading GT.

Published observations are exposed replay data, not a new blind experiment.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'experiments')]
from physics import Model, ReducedKernel
from validate import CountingBook
from run_quotient_benchmark import algebraic, fullj, qualified, compare
from run_branch_study import write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    if out.exists():
        parser.error('Use a new output directory; existing evidence is never overwritten.')
    original = json.loads((ROOT / 'experiments/QUOTIENT_HOLDOUT_CONTRACT.json').read_text())
    contract = json.loads((ROOT / 'experiments/QUOTIENT_HOLDOUT_RUN_CONTRACT.json').read_text())
    public = ROOT / 'experiments/quotient_holdout_public'
    manifest = json.loads((public / 'manifest.json').read_text())
    out.mkdir(parents=True)
    outputs, receipts = [], []
    for method in contract['methods']:
        begin = time.perf_counter()
        book = CountingBook()
        model = Model(book, **original['geometry'])
        kernel = ReducedKernel(model, np.full(4, .5), active=np.arange(4))
        for entry in manifest:
            case = entry['case']
            with np.load(public / entry['observations']) as data:
                y = data['y']
            row = dict(id=case, case=case, y=y, box=original['box'],
                       fixed=np.full(4, 15.25),
                       starts=np.asarray(contract['starts64'] if method == 'ReducedFullJ64'
                                         else contract['starts8']))
            result = (algebraic(model, kernel, row, method)
                      if method in ('quotient', 'six') else fullj(model, row, kernel))
            result.update(case=case, method=method, exposure='EXPOSED_REPLAY_NOT_NEW_BLIND')
            write(out / f'{case}_{method}.json', result)
            outputs.append(result)
            print(case, method, result['status'], 'qualified', len(qualified(result)), flush=True)
        receipt = dict(method=method, total_wall_including_target_logging=time.perf_counter()-begin,
                       accounting=book.summary(model.ns, model.nr))
        write(out / f'{method}_accounting.json', receipt)
        receipts.append(dict(method=method,
                             total_wall_seconds=receipt['total_wall_including_target_logging'],
                             full_complex_RHS=receipt['accounting']['total_complex_rhs']))
    cases = []
    for entry in manifest:
        selected = {r['method']: r for r in outputs if r['case'] == entry['case']}
        cases.append(dict(case=entry['case'],
                          qualified_counts={m: len(qualified(r)) for m, r in selected.items()},
                          quotient_vs_controls={m: compare(selected['quotient'], selected[m])
                                                for m in selected if m != 'quotient'}))
    write(out / 'SUMMARY.json', dict(status='EXPOSED_REPLAY_COMPLETE', truth_access='NONE',
          original_generation_not_repeated=True, methods=receipts, cases=cases,
          new_full_complex_RHS=sum(r['full_complex_RHS'] for r in receipts)))


if __name__ == '__main__':
    main()
