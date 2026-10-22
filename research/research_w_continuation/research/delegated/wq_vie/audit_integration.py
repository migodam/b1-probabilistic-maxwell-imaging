"""Non-solve audit: singular integration and the unchanged geometry sampler."""
from pathlib import Path
from types import SimpleNamespace
import json
import time

import numpy as np
from scipy.special import hankel1

from p0.geometry import material_grid, probes, target_centroid
from p0.vie import square_self_green
from wq.vie_upgrade import (
    rectangle_self_green, square_averaged_self_green,
    offset_source_integrals, frozen_subcell_material,
)


def grid(nx, frequency=1.):
    h = 4/nx
    ny = 3*nx//4
    x = -2+(np.arange(nx)+.5)*h
    y = (np.arange(ny)+.5)*h
    xx, yy = np.meshgrid(x, y, indexing='xy')
    return SimpleNamespace(nx=nx, ny=ny, N=nx*ny, h=h,
                           freq=frequency, points=np.c_[xx.ravel(), yy.ravel()])


def duffy_rectangle(k, h, offset, order=192):
    """Independent tensor quadrature after a two-triangle Duffy map."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    u = (nodes[:, None]+1)/2
    v = (nodes[None, :]+1)/2
    weights = weights[:, None]*weights[None, :]/4
    result = 0j
    for lx in (h/2-offset[0], h/2+offset[0]):
        for ly in (h/2-offset[1], h/2+offset[1]):
            if lx == 0 or ly == 0:
                continue
            for x, y in ((lx*u, ly*u*v), (lx*u*v, ly*u)):
                green = .25j*hankel1(0, k*np.hypot(x, y))
                result += np.sum(weights*lx*ly*u*green)
    return result


def duffy_double_square(k, h, order=192):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    u = (nodes[:, None]+1)/2
    v = (nodes[None, :]+1)/2
    weights = weights[:, None]*weights[None, :]/4
    # Two triangles per quadrant.  Symmetry gives eight identical terms.
    return 8*np.sum(weights*h*h*u*(1-u)*(1-u*v)
                    *.25j*hankel1(0, k*h*u*np.sqrt(1+v*v)))


def half_cell(k, h, destination, side, order=64):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x = side*h/4+nodes[:, None]*h/4
    y = nodes[None, :]*h/2
    radius = np.hypot(destination[0]-x, destination[1]-y)
    return np.sum(weights[:, None]*weights[None, :]*h*h/8
                  *.25j*hankel1(0, k*radius))


def finite_centroid(value):
    return None if not np.isfinite(value).all() else value.tolist()


def main():
    started = time.process_time()
    wall_started = time.perf_counter()
    out = Path(__file__).resolve().parent
    h = 4/96
    k = 2*np.pi
    checks = []
    corner = rectangle_self_green(k, h, (h/2, h/2))
    expected = square_self_green(k, 2*h)/4
    checks.append(dict(name='corner_quarter_of_double_square',
                       relative_error=float(abs(corner-expected)/abs(expected))))
    offset = np.array([-.2*h, .13*h])
    reference = duffy_rectangle(k, h, offset)
    value = rectangle_self_green(k, h, offset)
    checks.append(dict(name='offcenter_rectangle_vs_independent_duffy192',
                       relative_error=float(abs(value-reference)/abs(reference))))
    shifts = h*np.array([(a, b) for a in [-.375, -.125, .125, .375]
                         for b in [-.375, -.125, .125, .375]])
    partition = offset_source_integrals(-shifts, k, h/4, 16).sum()
    expected = square_self_green(k, h)
    checks.append(dict(name='four_by_four_microcell_partition_of_self',
                       relative_error=float(abs(partition-expected)/abs(expected))))
    averaged = square_averaged_self_green(k, h)
    reference = duffy_double_square(k, h)
    checks.append(dict(name='double_square_self_vs_independent_duffy192',
                       relative_error=float(abs(averaged-reference)/abs(reference)),
                       value_real=float(averaged.real), value_imag=float(averaged.imag)))
    destination = np.array([h, 0.])
    contrast = np.where(shifts[:, 0] < 0, 1+.01j, 3+.02j)
    weighted = np.dot(offset_source_integrals(destination-shifts, k, h/4, 8), contrast)
    reference = ((1+.01j)*half_cell(k, h, destination, -1)
                 + (3+.02j)*half_cell(k, h, destination, 1))
    mean_only = contrast.mean()*offset_source_integrals(destination[None], k, h, 32)[0]
    checks.append(dict(name='contrast_weighted_half_cell_vs_split_rectangle',
                       relative_error=float(abs(weighted-reference)/abs(reference)),
                       mean_chi_shortcut_relative_difference=float(abs(mean_only-reference)/abs(reference))))
    integration_cpu = time.process_time()-started

    rows = []
    geometry_started = time.process_time()
    geometry = grid(96)
    for seed in (2026102201, 2026102202):
        for label, probe in (('cal', None), ('task', probes()[8])):
            previous = None
            for sub in (4, 8, 16, 32):
                epsilon = material_grid(geometry, seed, subsamples=sub, probe=probe)
                chi = epsilon-1
                rows.append(dict(seed=seed, label=label, nx=96, subsamples=sub,
                                 active_cells=int(np.count_nonzero(chi)),
                                 dielectric_mass=float(chi.real.sum()*geometry.h**2),
                                 conductivity_mass=float(chi.imag.sum()*geometry.h**2),
                                 area_L2_epsilon_change=None if previous is None else
                                 float(np.linalg.norm(epsilon-previous)*geometry.h),
                                 max_epsilon_change=None if previous is None else
                                 float(np.max(np.abs(epsilon-previous))),
                                 target_centroid=finite_centroid(target_centroid(geometry.points, epsilon))))
                previous = epsilon
    geometry_cpu = time.process_time()-geometry_started
    sample_started = time.process_time()
    sample_grid = grid(32)
    sample_rows = []
    for sub in (4, 8, 16, 32):
        material = frozen_subcell_material(sample_grid, 2026102201, sub, probes()[8])
        old = material_grid(sample_grid, 2026102201, subsamples=sub, probe=probes()[8])
        sample_rows.append(dict(nx=32, subsamples=sub,
                                epsilon_array_exactly_equal_old=bool(np.array_equal(material.epsilon, old)),
                                max_contrast_mean_roundoff=float(np.max(np.abs(
                                    1+material.contrast.mean(axis=0)-old[material.active]))),
                                contrast_storage_bytes=int(material.contrast.nbytes)))
    sample_cpu = time.process_time()-sample_started
    summary = dict(status='INTEGRATION_AND_GEOMETRY_ONLY_NO_PHYSICAL_SOLVES',
                   physical_rhs=0, physical_lu=0,
                   integration_cpu_seconds=integration_cpu,
                   geometry_audit_cpu_seconds=geometry_cpu,
                   sampled_material_cpu_seconds=sample_cpu,
                   audit_cpu_seconds=time.process_time()-started,
                   total_process_cpu_seconds=time.process_time(),
                   wall_seconds=time.perf_counter()-wall_started,
                   integration_checks=checks,
                   material_identity_checks=sample_rows,
                   geometry_rows=rows,
                   limitations=['Duffy192 is a numerical integration reference, not a PDE certificate.',
                                'No field approximation, Maxwell convergence or qualification gate was tested.'])
    (out/'NON_SOLVE_AUDIT.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: summary[key] for key in
                      ('status', 'physical_rhs', 'physical_lu', 'total_process_cpu_seconds',
                       'wall_seconds', 'integration_checks', 'material_identity_checks')}, indent=2))


if __name__ == '__main__':
    main()
