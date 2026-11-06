"""One public-reference rank diagnostic; no truth, training, or new NN arm."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY.parents[1]))
from research.research_z_alg_r4_diversity.src.physics import Model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--anchor', type=Path)
    parser.add_argument('--output', type=Path, default=STUDY / 'audit/READOUT_GEOMETRY.json')
    args = parser.parse_args()
    started = time.perf_counter()
    with np.load(STUDY / 'data/operators.npz') as z:
        public = {k: z[k] for k in z.files}
    rows = []
    with threadpool_limits(limits=1):
        model = Model(public)
        if args.anchor:
            with np.load(args.anchor) as z:
                m, total = z['m0'], z['Q0']
        else:
            m = np.full(576, .1)
            total = model.forward(m, label='public_reference_readout_geometry')['Q']
        L = np.eye(576) - m[:, None] * model.G
        LV = L @ public['V48']
        LP0Lstar = L @ L.conj().T - LV @ LV.conj().T
        for p in (1, 3, 6, 12):
            indices = np.arange(0, 12, 12 // p)
            Q = total[:, indices]
            d = np.sum(np.abs(Q) ** 2, axis=1)
            C = np.real(LP0Lstar * (Q @ Q.conj().T).conj()) / d[:, None] / d[None, :]
            eigenvalues = eigvalsh(C)
            adjoint_field = L.conj().T @ Q
            null_adjoint = adjoint_field - public['V48'] @ (public['V48'].conj().T @ adjoint_field)
            rows.append(dict(sources=p, source_indices=indices.tolist(),
                             material_dimension=576,
                             null_current_real_dimension=2 * 528 * p,
                             gram_eigenvalues=eigenvalues.tolist(),
                             eigenvalue_min=float(eigenvalues[0]),
                             eigenvalue_max=float(eigenvalues[-1]),
                             eigenvalue_ratio=float(eigenvalues[0] / eigenvalues[-1]),
                             numerical_rank_relative_1e_10=int(np.sum(eigenvalues > eigenvalues[-1] * 1e-10)),
                             numerical_rank_relative_1e_6=int(np.sum(eigenvalues > eigenvalues[-1] * 1e-6)),
                             weighted_mean_annihilator_relative=float(np.linalg.norm(null_adjoint) / np.linalg.norm(adjoint_field)),
                             weighted_mean_direction='h_i=sum_s |Q_is|^2; conditional annihilator, not claimed exact at finite receiver sampling',
                             field_row_energy_min=float(d.min())))
    receipt = dict(status='PUBLIC_REFERENCE_DIAGNOSTIC', truth_access='NONE', reference_material=.1 if args.anchor is None else None,
                   anchor_input=str(args.anchor) if args.anchor else None,
                   wall_seconds=time.perf_counter() - started, counts=dict(model.counts),
                   rows=rows, interpretation='Readout differential rank is not Maxwell material observability or NN learnability; no threshold truncation enters training.')
    out = args.output
    out.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({**{k: v for k, v in receipt.items() if k != 'rows'},
                      'rows': [{k: v for k, v in r.items() if k != 'gram_eigenvalues'} for r in rows]}, indent=2))


if __name__ == '__main__':
    main()
