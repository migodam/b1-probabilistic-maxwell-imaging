"""Independent Fourier–Bessel replay of the 30-wavelength two-slot construction."""

import os
for name in ('OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[name] = '1'

from collections import Counter
from itertools import product
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import h1vp, hankel1, jv, jvp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from separated import encoded, l2, newton, pack, root_gaps, write_json

K, RADIUS, SPACING = 2*np.pi, .04, 30
ORDERS = (0, 1, 2, 3, 4, 6, 8, 12, 16)
BASE = np.array([2., 8.])


def coordinates(x, y):
    offset = x[:, None, :] - y[None, :, :]
    return np.linalg.norm(offset, axis=-1), np.angle(offset[..., 0] + 1j*offset[..., 1]), offset


def relative(a, b):
    return l2(a-b) / max(l2(b), 1e-300)


class Cylinders:
    """State unknowns are outgoing partial waves normalized at each circular interface."""

    def __init__(self, components, spacing, order, rx_x=.013, *, k=K, direction=None):
        if direction is None:
            shifts = np.c_[spacing*np.arange(components), np.zeros(components)]
        else:
            shifts = spacing*np.arange(components)[:, None]*np.asarray(direction)[None, :]
        self.centers = (shifts[:, None, :] + np.array([[-.18, 0.], [.18, 0.]])).reshape(-1, 2)
        self.sources = shifts + np.array([0., -.85])
        self.receivers = shifts + np.array([rx_x, .95])
        self.n = np.arange(-order, order+1)
        self.modes = len(self.n)
        self.disks = len(self.centers)
        self.spacing = spacing
        self.k = k
        self.identity = np.eye(self.disks*self.modes, dtype=complex)
        self.local = np.eye(components, dtype=bool)
        self.lo, self.hi = .5, 30.
        self.cost = Counter()
        self.target = np.zeros(components, dtype=complex)
        self.Jb = jv(self.n, self.k*RADIUS)
        self.Hb = hankel1(self.n, self.k*RADIUS)
        self.Jb_r = self.k*jvp(self.n, self.k*RADIUS)
        self.Hb_r = self.k*h1vp(self.n, self.k*RADIUS)

        dist, angle, _ = coordinates(self.centers, self.centers)
        np.fill_diagonal(dist, 1.)
        p = self.n[None, :] - self.n[:, None]
        G = (hankel1(p, self.k*dist[:, :, None, None]) * np.exp(1j*angle[:, :, None, None]*p)
             * self.Jb[:, None] / self.Hb[None, :])
        G[np.arange(self.disks), np.arange(self.disks)] = 0
        self.G = G.transpose(0, 2, 1, 3).reshape(self.identity.shape)

        dist, angle, _ = coordinates(self.centers, self.sources)
        B = .25j*hankel1(-self.n, self.k*dist[:, :, None])*np.exp(-1j*angle[:, :, None]*self.n)*self.Jb
        self.B = B.transpose(0, 2, 1).reshape(self.disks*self.modes, components)
        self.S = self.outgoing(self.receivers).reshape(components, self.disks*self.modes)
        self.cost['operator_assemblies'] += 1

    def outgoing(self, points):
        dist, angle, _ = coordinates(points, self.centers)
        return hankel1(self.n, self.k*dist[:, :, None])*np.exp(1j*angle[:, :, None]*self.n)/self.Hb

    def transmission(self, contrast, differentiate):
        scale = np.sqrt(1 + contrast[:, None])
        ki = self.k*scale
        z = ki*RADIUS
        ji, jbi = jv(self.n, z), jvp(self.n, z)
        ji_r = ki*jbi
        numerator = ji_r*self.Jb - ji*self.Jb_r
        denominator = ji*self.Hb_r - ji_r*self.Hb
        normalized = numerator/denominator * self.Hb/self.Jb
        if not differentiate:
            return normalized.ravel(), None
        d_ki = .5*self.k/scale
        d_ji = (RADIUS*d_ki)*jbi
        d_ji_r = d_ki*(jbi + z*jvp(self.n, z, n=2))
        dn = d_ji_r*self.Jb - d_ji*self.Jb_r
        dd = d_ji*self.Hb_r - d_ji_r*self.Hb
        derivative = (dn - numerator/denominator*dd)/denominator*self.Hb/self.Jb
        return normalized.ravel(), derivative.ravel()

    def solve(self, contrast, differentiate=False):
        tau, d_tau = self.transmission(contrast, differentiate)
        A = self.identity - tau[:, None]*self.G
        b = tau[:, None]*self.B
        factored = lu_factor(A, check_finite=False)
        a = lu_solve(factored, b, check_finite=False)
        self.cost.update({'state_LU': 1, 'source_RHS': b.shape[1]})
        y = self.S @ a
        f = pack(y[self.local] - self.target)
        if not differentiate:
            return y, f
        incoming = self.B + self.G @ a
        receiver = lu_solve(factored, self.S.T, trans=1, check_finite=False).T
        self.cost['receiver_transpose_RHS'] += self.S.shape[0]
        port_derivative = receiver[:, None, :]*(d_tau[:, None]*incoming).T[None, :, :]
        J = pack(port_derivative[self.local].reshape(len(self.target), self.disks, self.modes).sum(axis=2))
        return y, f, J, A, b, a

    def explicit_boundary_check(self, material, samples=48):
        """Evaluate exterior fields with direct Hankel functions and test both jumps."""
        _, _, _, _, _, coefficients = self.solve(material, differentiate=True)
        incoming = self.B + self.G @ coefficients
        modes = np.exp(1j*np.linspace(0., 2*np.pi, samples, endpoint=False)[:, None]*self.n)
        theta = np.linspace(0., 2*np.pi, samples, endpoint=False)
        normal = np.c_[np.cos(theta), np.sin(theta)]
        interior_at_boundary = (incoming + coefficients).reshape(self.disks, self.modes, -1)
        ki = self.k*np.sqrt(1 + material[:, None])
        derivative_ratio = ki*jvp(self.n, RADIUS*ki)/jv(self.n, RADIUS*ki)
        jumps = []
        for idx, center in enumerate(self.centers):
            points = center + RADIUS*normal
            src_dist, _, src_diff = coordinates(points, self.sources)
            src_proj = np.einsum('pki,pi->pk', src_diff/src_dist[:, :, None], normal)
            direct_incident = .25j*hankel1(0, self.k*src_dist)
            direct_incident_radial = .25j*self.k*h1vp(0, self.k*src_dist)*src_proj
            distance, angle, offset = coordinates(points, self.centers)
            radial = offset / distance[:, :, None]
            ortho = np.stack((-radial[:, :, 1], radial[:, :, 0]), axis=2)
            radial_proj = np.einsum('pdi,pi->pd', radial, normal)
            angular_proj = np.einsum('pdi,pi->pd', ortho, normal) / distance
            phase = np.exp(1j*angle[:, :, None]*self.n)
            wave = hankel1(self.n, self.k*distance[:, :, None])
            dw = (self.k*h1vp(self.n, self.k*distance[:, :, None])*radial_proj[:, :, None]
                  + 1j*self.n*wave*angular_proj[:, :, None])
            exterior = direct_incident + (wave*phase/self.Hb).reshape(samples, -1) @ coefficients
            exterior_derivative = direct_incident_radial + (dw*phase/self.Hb).reshape(samples, -1) @ coefficients
            inner = interior_at_boundary[idx]
            interior = modes @ inner
            interior_derivative = modes @ (derivative_ratio[idx, :, None]*inner)
            jumps.append([relative(exterior, interior), relative(exterior_derivative, interior_derivative)])
        return dict(max_relative_field_jump=max(row[0] for row in jumps),
                    max_relative_normal_derivative_jump=max(row[1] for row in jumps), angles=samples)


def inspect(order, pulse, saved):
    start, cpu = perf_counter(), process_time()
    one = Cylinders(1, 0, order)
    y_a, _ = one.solve(BASE)
    one.target = y_a[one.local]
    b, _, _, branch_b = newton(one, BASE[::-1])
    if branch_b['status'] == 'numerical_root' and l2(b-BASE) > 1e-5:
        seeds = {'A': BASE, 'B': b}
    else:
        seeds = {}
    mirror = Cylinders(1, 0, order, rx_x=0.)
    mirror_error = l2(mirror.solve(BASE)[0] - mirror.solve(BASE[::-1])[0])
    info = {'angular_order': order, 'retained_modes_per_disk': 2*order+1,
            'target': encoded(y_a), 'mirror_absolute_port_error': mirror_error,
            'local_partner': branch_b, 'local_partner_distinct': bool(seeds),
            'difference_from_pulse_local_target': l2(y_a-pulse['target']),
            'difference_from_pulse_local_partner_l2': l2(b-pulse['partner']), 'coupled_cases': []}
    for count in (2, 3):
        model = Cylinders(count, SPACING, order)
        model.target = np.repeat(one.target, count)
        roots, attempts = [], {}
        if seeds:
            for tag in product('AB', repeat=count):
                branch = ''.join(tag)
                seed = np.concatenate([seeds[k] for k in tag])
                z, y, J, attempt = newton(model, seed)
                attempts[branch] = attempt
                if attempt['status'] == 'numerical_root' and all(l2(z-r[1]) > 1e-5 for r in roots):
                    roots.append((branch, z, y))
                    saved[f'M{order}_L{count}_{branch}__m'] = z
                    saved[f'M{order}_L{count}_{branch}__ports'] = y
                    saved[f'M{order}_L{count}_{branch}__J'] = J
        gaps = root_gaps(roots, model.local)
        dx = [l2(z-pulse[count][tag]) for tag, z, _ in roots]
        boundary = model.explicit_boundary_check(roots[0][1]) if roots else None
        info['coupled_cases'].append({'components': count, 'numerical_distinct_roots': len(roots),
                                     'attempts': attempts, 'pairwise_gaps': gaps,
                                     'max_root_displacement_from_matching_pulse_branch': max(dx, default=None),
                                     'explicit_boundary_check': boundary, 'cost': dict(model.cost)})
    info.update(wall_seconds=perf_counter()-start, cpu_seconds=process_time()-cpu,
                single_component_cost=dict(one.cost), mirror_cost=dict(mirror.cost))
    return info


def main():
    beginning, cpu = perf_counter(), process_time()
    data = json.loads((HERE.parent / 'conditioning/results.json').read_text(encoding='utf-8'))
    angles = data['displaced_pair']
    re, im = angles['A']['physical_scattered_port'][0]
    pulse = {'target': re+1j*im, 'partner': np.array(angles['B']['material'])}
    for case in data['coupled_cases']:
        if case['spacing_reference_wavelengths'] == SPACING:
            pulse[case['components']] = {row['branch_seed']: np.array(row['material_at_stop'])
                                         for row in case['attempts']}
    saved = {}
    results = [inspect(order, pulse, saved) for order in ORDERS]
    for before, after in zip(results[:-1], results[1:], strict=True):
        prev_order, order = before['angular_order'], after['angular_order']
        for row in after['coupled_cases']:
            count = row['components']
            diffs = [l2(saved[f'M{order}_L{count}_{tag}__m']-saved[f'M{prev_order}_L{count}_{tag}__m'])
                     for tag in row['attempts'] if f'M{order}_L{count}_{tag}__m' in saved
                     and f'M{prev_order}_L{count}_{tag}__m' in saved]
            row['max_material_L2_change_from_previous_order'] = max(diffs, default=None)
            row['matched_roots_across_orders'] = len(diffs)
    np.savez_compressed(HERE / 'materials_and_ports.npz', **saved)
    metadata = {'scope': 'Fourier–Bessel truncation, lossless circular interfaces, no rigorous tail bound',
                'geometry_wavelength_units': {'disk_radius': RADIUS, 'centers_per_component': [[-.18, 0], [.18, 0]],
                    'source': [0, -.85], 'receiver': [.013, .95], 'component_spacing': SPACING},
                'fixed_trial_orders': list(ORDERS), 'runs': results,
                'wall_seconds': perf_counter()-beginning, 'cpu_seconds': process_time()-cpu}
    write_json(HERE / 'results.json', metadata)
    for result in results:
        print(result['angular_order'], result['local_partner']['status'],
              [(r['components'], r['numerical_distinct_roots']) for r in result['coupled_cases']], flush=True)


if __name__ == '__main__':
    main()
