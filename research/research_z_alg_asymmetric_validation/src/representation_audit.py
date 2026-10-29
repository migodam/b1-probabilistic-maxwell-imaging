"""E3 saved-material representation audit. No forward solver is imported or called.

Default invocation only describes readiness. Parent authorizes --execute after E1
qualified ENDPOINT.npz files exist. All representations use the SAME material
prior through change of variables, never independent normal coefficient priors.
"""
import argparse
import csv
import json
import os
from pathlib import Path

for _name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS']:
    os.environ[_name] = '1'

import numpy as np
from accounting import Book, P, write


def transforms(n=6):
    """Return real orthogonal coordinate maps, full complex DFT, and Haar tree."""
    if n != 6:
        raise ValueError('This frozen E3 implements exactly six real cells')
    index = np.arange(n)
    dft = np.exp(-2j * np.pi * np.outer(index, index) / n) / np.sqrt(n)
    # Independent REAL coordinates on the conjugate-symmetric complex DFT image.
    real_dft = np.stack([dft[0].real, np.sqrt(2) * dft[1].real,
                         np.sqrt(2) * dft[1].imag, np.sqrt(2) * dft[2].real,
                         np.sqrt(2) * dft[2].imag, dft[3].real])
    rows = [np.ones(n) / np.sqrt(n)]
    tree = []

    def split(lo, hi):
        if hi - lo < 2:
            return
        mid = (lo + hi) // 2
        left, right, total = mid - lo, hi - mid, hi - lo
        contrast = np.zeros(n)
        contrast[lo:mid] = np.sqrt(right / (left * total))
        contrast[mid:hi] = -np.sqrt(left / (right * total))
        tree.append(dict(row=len(rows), parent=[lo, hi], left=[lo, mid],
                         right=[mid, hi], interval_convention='half-open'))
        rows.append(contrast)
        split(lo, mid)
        split(mid, hi)

    split(0, n)
    return {'native': np.eye(n), 'DFT_real_conjugacy': real_dft,
            'unbalanced_Haar6': np.stack(rows)}, dft, tree


def load_materials():
    """Only qualified saved E1 endpoints; inspect result metadata, never truth."""
    records, materials = [], []
    for summary_path in sorted((P / 'results/E1').glob('P*/[AB]/SUMMARY.json')):
        summary = json.loads(summary_path.read_text())
        for row in summary['rows']:
            if not row.get('qualified', False) or not row.get('legal', False):
                continue
            source = summary_path.parent / row['key'] / 'ENDPOINT.npz'
            if not source.exists():
                raise RuntimeError('Qualified saved endpoint missing: ' + str(source))
            with np.load(source, allow_pickle=False) as saved:
                m = np.asarray(saved['m'])
                if np.iscomplexobj(m) and np.any(m.imag != 0):
                    raise ValueError('E3 requires real material: ' + str(source))
                m = np.asarray(m.real, dtype=np.float64)
            if m.shape != (6,) or not np.all(np.isfinite(m)):
                raise ValueError('Invalid six-cell saved material: ' + str(source))
            if np.any(m < .5) or np.any(m > 4.):
                raise ValueError('Claimed legal material outside frozen box')
            records.append(dict(index=len(records), parent=summary['parent'],
                                method=summary['method'], key=row['key'],
                                source=source.relative_to(P).as_posix(),
                                source_summary=summary_path.relative_to(P).as_posix()))
            materials.append(m)
    if not materials:
        raise RuntimeError('No completed E1 qualified saved six-cell endpoints')
    # Preserve repeated material/source occurrences; this is not posterior sampling.
    return np.stack(materials), records


def log_material_density(materials, lower=.5, upper=4., support_tolerance=0.):
    materials = np.asarray(materials)
    support = np.all((materials >= lower - support_tolerance)
                     & (materials <= upper + support_tolerance), axis=-1)
    density = np.full(materials.shape[:-1], -np.inf)
    density[support] = np.sum(np.log(2 * materials[support]
                                     / (upper ** 2 - lower ** 2)), axis=-1)
    return density, support


def execute():
    contract = json.loads((P / 'EXPERIMENT_CONTRACT.json').read_text())
    if contract['E1']['n_real_material'] != 6:
        raise ValueError('E3 frozen six-cell basis requires E1 n=6')
    if contract['new_prior']['name'] != 'P_Z_real_product_linear':
        raise ValueError('Unexpected prior; no silent prior replacement')
    # Check saved input availability before starting a single-use Book run.
    materials, sources = load_materials()
    out = P / 'results/E3'
    if out.exists():
        raise FileExistsError('E3 results already exist; preserve original run')
    book = Book('E3_saved_material_representation')
    actual_before, formula_before = book.actual, book.formula
    out.mkdir(parents=True, exist_ok=False)
    try:
        with book.op('E3_saved_material_basis_arithmetic', rhs=0, formula=0,
                     material_source_count=len(materials), forward_calls=0):
            maps, full_dft, tree = transforms()
            density, support = log_material_density(materials)
            arrays = dict(materials=materials, full_unitary_DFT=full_dft,
                          full_DFT_coefficients=materials @ full_dft.T,
                          native_log_prior=density)
            transform_report, rows = {}, []
            tol = 2e-12
            for name, transform in maps.items():
                coef = materials @ transform.T
                reconstructed = coef @ transform
                error = reconstructed - materials
                sign, logabsdet = np.linalg.slogdet(transform)
                restored_log_density, pushed_support = log_material_density(
                    reconstructed, support_tolerance=tol)
                # q(a) = p(U.T a) / |det U|, with the original box-image support.
                pushed_log_density = restored_log_density - logabsdet
                prior_log_error = np.abs(pushed_log_density - density)
                metric_error = np.abs(np.sum(coef ** 2, axis=1)
                                      - np.sum(materials ** 2, axis=1))
                gram_error = float(np.max(np.abs(transform @ transform.T - np.eye(6))))
                # A rank-three diagnostic: zeros in remaining independent REAL
                # coordinates. It is a projection, not an invertible new decoder.
                truncated_coef = coef.copy()
                truncated_coef[:, 3:] = 0.
                truncated = truncated_coef @ transform
                truncation_delta = truncated - materials
                _, truncated_legal = log_material_density(truncated)
                arrays.update({name + '_matrix': transform,
                               name + '_coefficients': coef,
                               name + '_roundtrip': reconstructed,
                               name + '_roundtrip_error': error,
                               name + '_log_prior_pushforward': pushed_log_density,
                               name + '_truncated3_material': truncated,
                               name + '_truncated3_error': truncation_delta})
                transform_report[name] = dict(
                    matrix=transform, determinant_sign=float(sign),
                    absolute_determinant=float(np.exp(logabsdet)),
                    log_absolute_determinant=float(logabsdet),
                    real_coordinate_dimension=6, orthogonality_max_abs=gram_error,
                    max_roundtrip_abs=float(np.max(np.abs(error))),
                    max_metric_squared_norm_error=float(np.max(metric_error)),
                    max_log_prior_density_error=float(np.max(prior_log_error)),
                    all_pushforward_support=bool(np.all(pushed_support)),
                    support_definition='a such that every component of U.T@a is in [0.5,4]',
                    pushforward_density='product_i(2*(U.T@a)_i/(4^2-0.5^2))/abs(det U) on U(box)',
                    independent_normal_coefficient_prior=False,
                    determinant_scope='Real six-dimensional induced coordinates; not ambient C6 volume',
                    floating_support_tolerance=tol,
                    truncated3=dict(rank=3, retained_independent_real_coordinates=[0, 1, 2],
                                    discarded_coordinates=[3, 4, 5],
                                    max_material_RMS=float(np.max(np.sqrt(np.mean(truncation_delta**2, axis=1)))),
                                    legal_projected_count=int(np.sum(truncated_legal)),
                                    support_changed=True,
                                    original_absolutely_continuous_prior_mass_of_exact_rank3_support=0.,
                                    physical_fiber_or_gauge_claim=False,
                                    note='No forward is computed; projection changes material and prior support even when projected material remains in box'))
                for i, source in enumerate(sources):
                    rows.append(dict(**source, transform=name,
                                     roundtrip_max_abs=float(np.max(np.abs(error[i]))),
                                     roundtrip_material_RMS=float(np.sqrt(np.mean(error[i]**2))),
                                     metric_squared_norm_error=float(metric_error[i]),
                                     native_log_prior=float(density[i]),
                                     pushed_log_prior=float(pushed_log_density[i]),
                                     log_prior_error=float(prior_log_error[i]),
                                     pushforward_support=bool(pushed_support[i]),
                                     truncated3_material_RMS=float(np.sqrt(np.mean(truncation_delta[i]**2))),
                                     truncated3_in_box=bool(truncated_legal[i])))

            fourier = arrays['full_DFT_coefficients']
            conjugacy_error = float(np.max(np.abs(fourier[:, [5, 4]]
                                                  - fourier[:, [1, 2]].conj())))
            real_mode_imag_error = float(np.max(np.abs(fourier[:, [0, 3]].imag)))
            complex_roundtrip = fourier @ full_dft.conj()
            complex_roundtrip_error = float(np.max(np.abs(complex_roundtrip - materials)))
            unitary_error = float(np.max(np.abs(full_dft @ full_dft.conj().T - np.eye(6))))
            gate = all(r['max_roundtrip_abs'] <= tol
                       and r['orthogonality_max_abs'] <= tol
                       and r['max_log_prior_density_error'] <= tol
                       and r['all_pushforward_support']
                       for r in transform_report.values())
            gate = bool(gate and conjugacy_error <= tol and real_mode_imag_error <= tol
                        and unitary_error <= tol and complex_roundtrip_error <= tol)
            np.savez_compressed(out / 'COMPACT_ARRAYS.npz', **arrays)
            with (out / 'PER_MATERIAL.csv').open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            summary = dict(
                status='PASS' if gate else 'FAIL',
                qualification='Complete representation roundtrip/prior arithmetic only',
                source_material_occurrences=len(materials), source_records=sources,
                input_scope='Qualified ENDPOINT.npz from completed E1 summaries only; repetitions preserved; not probability samples',
                new_forward_calls=0, new_actual_complex_RHS=0,
                new_closed_form_source_RHS_equivalent=0, transforms=transform_report,
                unbalanced_Haar_tree=tree,
                fourier=dict(unitary_max_abs=unitary_error,
                             conjugacy_max_abs=conjugacy_error,
                             zero_and_nyquist_imag_max_abs=real_mode_imag_error,
                             full_complex_roundtrip_max_abs=complex_roundtrip_error,
                             real_coordinates=['c0.real', 'sqrt(2)*c1.real', 'sqrt(2)*c1.imag',
                                               'sqrt(2)*c2.real', 'sqrt(2)*c2.imag', 'c3.real'],
                             ambient_complex_support='c5=conj(c1),c4=conj(c2);c0,c3 real',
                             complex_ambient_density='Singular six-real-dimensional support in twelve-real-dimensional C6; no ambient Lebesgue density'),
                same_prior='P_Z_real_product_linear pushed through each complete invertible map',
                same_material_metric='Euclidean norm/sqrt(6); real complete maps orthogonal and full complex DFT unitary',
                changed_prior_or_support_for_complete_maps=False,
                Full66_firewall='Original P1 Gaussian66 decoder prior is a different support/model; none of these bases silently covers or replaces it',
                truncation_scope='Three coefficients are a rank-three support-change diagnostic; no forward response, fiber, gauge, or imaging claim',
                exactness='Algebraic orthogonality/unitarity identities with floating roundtrip diagnostics; no interval arithmetic certificate')
            write(out / 'SUMMARY.json', summary)
        if (book.actual - actual_before, book.formula - formula_before) != (0, 0):
            raise RuntimeError('E3 unexpectedly incurred physics cost')
        book.end('COMPLETE' if gate else 'FAILED_REPRESENTATION_QA',
                 new_physics_RHS=0, output='results/E3/SUMMARY.json')
        print(json.dumps(dict(status=summary['status'], source_material_occurrences=len(materials),
                              new_physics_RHS=0, output=str(out))))
    except BaseException as exc:
        write(out / 'INTERRUPTED.json', dict(error=repr(exc), new_physics_RHS=0))
        book.end('FAILED', error=repr(exc), new_physics_RHS=0)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if not args.execute:
        print(json.dumps(dict(status='READY_NOT_EXECUTED', new_physics_RHS=0,
                              requires='Parent authorization and completed E1 qualified saved ENDPOINT.npz',
                              command='python src/representation_audit.py --execute')))
        return
    execute()


if __name__ == '__main__':
    main()
