"""Zero-Maxwell geometry audit. No solver, FEM, mesh constructor or LU imports.

The historical scalar ``minimum_feature`` does not define a geometric class.
``raster_disk_coverage`` is an explicitly NEW diagnostic contract: both phases
must be unions of closed Euclidean disks of radius minimum_feature/2, with
material extended by air outside the image. It is a raster diagnostic, not a
continuous certificate and not a retroactive definition of the frozen prior.
Exact sharp corners fail this two-sided rolling-disk contract.
"""
from __future__ import annotations
import argparse
import types
import json
from pathlib import Path
import time
import numpy as np
from scipy import ndimage
from scipy.optimize import brentq, differential_evolution

SEEDS = (2026102201, 2026102202)
DOMAIN = (-2., 2., 0., 3.)
QUERY = (-1.2, 1.2, 1.6, 2.8)


def phases(seed):
    return np.random.default_rng(seed).uniform(0., 2*np.pi, 5)


def wall_implicit(points, seed):
    """Same support equations as frozen wall_values, no coefficient smoothing."""
    x, y = np.asarray(points, dtype=float).T
    p = phases(seed)
    low = .64+.09*x+.055*np.sin(2.4*x+p[0])
    thick = .48+.085*np.cos(2.1*x+p[1])
    angle = np.arctan2(y-.92, x+.20)
    rad = np.hypot((x+.20)/.25, (y-.92)/.15)
    return np.minimum.reduce([1.48+.08*np.sin(4*y+p[2])-np.abs(x),
                              y-low, low+thick-y,
                              rad-1-.18*np.sin(3*angle+p[3])])


def cavity_boundary(theta, seed):
    """Exact polar parametrization, including the original physical atan2."""
    theta = np.asarray(theta)
    rho = ((1+.18*np.sin(3*theta+phases(seed)[3]))
           / np.sqrt((np.cos(theta)/.25)**2+(np.sin(theta)/.15)**2))
    return np.stack((-.20+rho*np.cos(theta), .92+rho*np.sin(theta)), axis=-1)


def strip_boundary(x, seed, side):
    p = phases(seed)
    low = .64+.09*x+.055*np.sin(2.4*x+p[0])
    return low if side == 'bottom' else low+.48+.085*np.cos(2.1*x+p[1])


def boundary_audit(seed):
    theta = np.linspace(0., 2*np.pi, 65537)
    c = cavity_boundary(theta, seed)
    result = {}
    for side in ('bottom', 'top'):
        def gap(t):
            z = cavity_boundary(t, seed)
            v = z[..., 1]-strip_boundary(z[..., 0], seed, side)
            return v if side == 'bottom' else -v
        g = gap(theta)
        roots = [float(brentq(gap, theta[i], theta[i+1], xtol=1e-14))
                 for i in np.flatnonzero(g[:-1]*g[1:] < 0)]
        minimum = differential_evolution(lambda a: float(gap(a[0])),
                                         [(0., 2*np.pi)], seed=0,
                                         tol=1e-11, polish=True)
        def objective(v):
            z = cavity_boundary(v[0], seed)
            q = np.array([v[1], strip_boundary(v[1], seed, side)])
            return float(np.sum((z-q)**2))
        nearest = differential_evolution(objective,
                    [(0., 2*np.pi), (-.6, .2)], seed=0,
                    tol=1e-11, polish=True)
        z = cavity_boundary(nearest.x[0], seed)
        q = np.array([nearest.x[1], strip_boundary(nearest.x[1], seed, side)])
        crossing_angles = []
        for t in roots:
            z = cavity_boundary(t, seed)
            dt = 1e-6
            tangent = (cavity_boundary(t+dt,seed)-cavity_boundary(t-dt,seed))/(2*dt)
            p = phases(seed)
            slope = .09+.132*np.cos(2.4*z[0]+p[0])
            if side == 'top':
                slope -= .1785*np.sin(2.1*z[0]+p[1])
            sine = abs(tangent[1]-slope*tangent[0])/(np.linalg.norm(tangent)*np.hypot(1,slope))
            crossing_angles.append(float(np.degrees(np.arcsin(np.clip(sine,0,1)))))
        result[side] = dict(signed_vertical_gap_min=float(minimum.fun),
            transverse_crossing_acute_angles_degrees=crossing_angles,
            intersection_theta=roots,
            intersection_points=[cavity_boundary(t, seed).tolist() for t in roots],
            closest_euclidean_gap=float(np.sqrt(max(0., nearest.fun))),
            closest_cavity_point=z.tolist(), closest_strip_point=q.tolist(),
            distance_status='numerical global optimizer diagnostic, not certified bound')
    result['cavity_topology'] = ('open_notch' if result['top']['intersection_theta']
                                  else 'enclosed_hole')
    result['thin_wedge_observation'] = (
        'At transverse cavity/strip intersection material wedges have arbitrarily '
        'small local widths; two-sided rolling-disk contract fails there.'
        if result['cavity_topology'] == 'open_notch' else None)
    return result


def contour_area(seed, resolution=1):
    """Decoder-only contour, no polygon mesher or PDE import."""
    import contourpy
    x = np.linspace(-1.65, 1.65, 1320*resolution+1)
    y = np.linspace(.3, 1.45, 576*resolution+1)
    xx, yy = np.meshgrid(x, y)
    phi = wall_implicit(np.c_[xx.ravel(), yy.ravel()], seed).reshape(xx.shape)
    loops = contourpy.contour_generator(x=x, y=y, z=phi).lines(0.)
    areas = []
    for a in loops:
        if np.linalg.norm(a[0]-a[-1]) > 1e-9:
            raise ValueError('Open contour in finite contour window')
        areas.append(float(.5*np.sum(a[:-1, 0]*a[1:, 1]-a[1:, 0]*a[:-1, 1])))
    # contourpy orients level sets consistently; largest outer contour and holes
    # have opposite signs. Absolute signed sum supports both orientations.
    return dict(resolution=resolution, dx=float(x[1]-x[0]),
                dy=float(y[1]-y[0]), loops=len(loops),
                occupied_area=abs(sum(areas)), signed_loop_areas=areas,
                status='polygon approximation; convergence diagnostic only')


def raster_disk_coverage(mask, spacing, minimum_feature=.1, *, pad_width=None):
    """Check NEW two-phase disk coverage on cell-center binary samples.

    spacing is (dy,dx). Material is false outside the supplied rectangle.
    A phase is opened by erosion and dilation with radius r=f/2 via Euclidean
    center distances. Pixels removed by opening are reported, never repaired.
    Center-to-cell uncertainty is half the cell diagonal; EDT reports distances
    between sample centers, not exact continuous-boundary distances. Nonzero
    losses prove failure of this discrete contract only. Even zero losses cannot
    certify the unsampled continuous geometry. No area tolerance is hidden.
    """
    mask = np.asarray(mask)
    if mask.ndim != 2 or mask.dtype != bool:
        raise ValueError('mask must be a two-dimensional boolean array')
    spacing = np.asarray(spacing, dtype=float)
    if spacing.shape != (2,) or np.any(spacing <= 0) or not np.isfinite(spacing).all():
        raise ValueError('spacing must contain finite positive dy,dx')
    if not np.isfinite(minimum_feature) or minimum_feature <= 0:
        raise ValueError('minimum_feature must be positive and finite')
    r = float(minimum_feature)/2
    n = int(np.ceil(3*r/spacing.min()))+4 if pad_width is None else int(pad_width)
    if n*spacing.min() < 2*r:
        raise ValueError('padding must span at least two radii')
    padded = np.pad(mask, n, constant_values=False)
    losses = {}
    for name, phase in (('material', padded), ('air', ~padded)):
        core = ndimage.distance_transform_edt(phase, sampling=spacing) > r
        if core.any():
            opened = ndimage.distance_transform_edt(~core, sampling=spacing) <= r
        else:
            opened = np.zeros_like(core)
        lost = (phase & ~opened)[n:-n, n:-n]
        losses[name] = dict(removed_pixels=int(lost.sum()),
                           removed_area=float(lost.sum()*np.prod(spacing)))
    return dict(contract='NEW_TWO_PHASE_EUCLIDEAN_DISK_COVERAGE_RASTER_V1',
                minimum_feature=float(minimum_feature), disk_radius=r,
                spacing_dy_dx=spacing.tolist(),
                sample_to_cell_radius=float(np.linalg.norm(spacing)/2),
                original_prior_qualification='UNDETERMINED_UNSPECIFIED_CRITERION',
                discrete_contract_pass=all(v['removed_pixels'] == 0 for v in losses.values()),
                losses=losses,
                continuous_certificate=False,
                notes='Sharp corners fail this stronger contract; no reinterpretation of frozen minimum_feature.')


def erosion_topology(mask, spacing, minimum_feature=.1, *, connectivity=8):
    """NEW corner-permitting topology-persistence diagnostic, no smoothing.

    Both phases use explicitly selected 4- or 8-neighbor connectivity. Every connected
    component must contain exactly one connected component after Euclidean
    erosion by r=f/2: no extinction or splitting. Air extends outside D.
    This controls topology-carrying necks/components, not every local width:
    thin spurs or redundant narrow gaps can survive this diagnostic. It is
    deliberately weaker than two-sided opening equality, and not the original
    scalar prior. Return three erosion radii to expose center-distance risk.
    """
    mask = np.asarray(mask)
    spacing = np.asarray(spacing,dtype=float)
    if mask.ndim != 2 or mask.dtype != bool or spacing.shape != (2,):
        raise ValueError('Expected boolean 2D mask and dy,dx')
    if not np.isfinite(spacing).all() or np.any(spacing <= 0) or not np.isfinite(minimum_feature) or minimum_feature <= 0:
        raise ValueError('Finite positive spacing and feature size required')
    r = float(minimum_feature)/2
    eta = float(np.linalg.norm(spacing)/2)
    pad = int(np.ceil(3*r/spacing.min()))+4
    padded = np.pad(mask,pad,constant_values=False)
    if connectivity not in (4,8):
        raise ValueError('connectivity must be 4 or 8')
    face = ndimage.generate_binary_structure(2,1 if connectivity == 4 else 2)
    result = {}
    for name,phase in (('material',padded),('air',~padded)):
        labels,n = ndimage.label(phase,structure=face)
        distances = ndimage.distance_transform_edt(phase,sampling=spacing)
        radii=[]
        for radius in (max(0.,r-eta),r,r+eta):
            eroded = distances > radius
            reduced,m = ndimage.label(eroded,structure=face)
            slices = ndimage.find_objects(reduced)
            inherited = np.zeros(m,dtype=int)
            for j,sl in enumerate(slices):
                loc = np.argwhere(reduced[sl] == j+1)[0]
                inherited[j] = labels[tuple(loc[k]+sl[k].start for k in range(2))]
            counts = np.bincount(inherited,minlength=n+1)[1:]
            radii.append(dict(radius=float(radius),original_components=int(n),
                eroded_components=int(m),extinct_components=int(np.count_nonzero(counts==0)),
                split_components=int(np.count_nonzero(counts>1)),
                pass_components=bool(np.all(counts==1))))
        result[name]=radii
    nominal=all(v[1]['pass_components'] for v in result.values())
    stable=all(all(z['pass_components']==v[1]['pass_components'] for z in v) for v in result.values())
    return dict(contract='PROSPECTIVE_W_R_GEOM_V1_TOPOLOGY_ONLY',
        minimum_feature=float(minimum_feature),disk_radius=r,connectivity=connectivity,
        spacing_dy_dx=spacing.tolist(),radius_uncertainty_probe=eta,
        nominal_discrete_pass=bool(nominal),radius_probe_stable=bool(stable),
        phase_components=result,continuous_certificate=False,
        original_prior_qualification='UNRESOLVED_SCALAR_SEMANTICS',
        limitations='Corners allowed; thin spurs/redundant narrow gaps are NOT excluded. No sampled topology result certifies continuous minimum width.')


def validate_material_samples(points, epsilon_inf, conductivity, *, query=QUERY):
    """Pointwise bounds and forced-query-air check; no unsampled certificate."""
    z = np.asarray(points, dtype=float)
    e, s = np.asarray(epsilon_inf), np.asarray(conductivity)
    if z.ndim != 2 or z.shape[1] != 2 or e.shape != (len(z),) or s.shape != (len(z),):
        raise ValueError('Expected points(n,2), epsilon_inf(n), conductivity(n)')
    if not (np.isfinite(z).all() and np.isfinite(e).all() and np.isfinite(s).all()):
        raise ValueError('Nonfinite input')
    ql, qr, qb, qt = query
    q = (z[:, 0] >= ql) & (z[:, 0] <= qr) & (z[:, 1] >= qb) & (z[:, 1] <= qt)
    occupied = (e != 1) | (s != 0)
    outside = ((z[:, 0] < DOMAIN[0]) | (z[:, 0] > DOMAIN[1])
               | (z[:, 1] < DOMAIN[2]) | (z[:, 1] > DOMAIN[3]))
    return dict(sample_count=len(z), bounds_pass=bool(((e >= 1)&(e <= 9)&(s >= 0)&(s <= .15)).all()),
                query_air_pass=bool((~occupied[q]).all()),
                outside_domain_air_pass=bool((~occupied[outside]).all()),
                scope='sampled points only; no geometry or minimum-feature certificate')


def disk_in_rectangle(center, radius, box=QUERY, *, strict=False):
    """Exact axis-aligned disk containment. strict=True adds NEW no-contact rule."""
    c = np.asarray(center, dtype=float)
    if c.shape != (2,) or not np.isfinite(c).all() or not np.isfinite(radius) or radius <= 0:
        raise ValueError('Finite center(2) and positive radius required')
    box = np.asarray(box, dtype=float)
    if box.shape != (4,) or not np.isfinite(box).all() or box[0] >= box[1] or box[2] >= box[3]:
        raise ValueError('box must be finite ordered xl,xr,yb,yt')
    margins = np.array([c[0]-radius-box[0], box[1]-c[0]-radius,
                        c[1]-radius-box[2], box[3]-c[1]-radius])
    return dict(contained=bool((margins > 0).all() if strict else (margins >= 0).all()),
                minimum_margin=float(margins.min()), strict_clearance=bool(strict))


def adversarial_cases():
    h = .0025
    x = np.arange(-.5+h/2, .5, h)
    xx, yy = np.meshgrid(x, x)
    body = (abs(xx) < .3)&(abs(yy) < .3)
    cases = {
        'thin_island_width_0.06': (abs(xx) < .03)&(abs(yy) < .2),
        'thin_air_channel_width_0.06': body & ~((((xx+.14)**2+yy**2 < .10**2)|((xx-.14)**2+yy**2 < .10**2)|((abs(xx)<.14)&(abs(yy)<.03)))),
        'redundant_thin_air_channel_counterexample': body & ~((abs(xx) < .03)&(yy > -.31)),
        'thin_hole_diameter_0.06': body & ~(xx**2+yy**2 < .03**2),
        'sharp_square_corners': body,
        'round_disk_diameter_0.24': xx**2+yy**2 < .12**2,
        'thin_neck_material_dumbbell': (((xx+.18)**2+yy**2 < .13**2)|((xx-.18)**2+yy**2 < .13**2)|((abs(xx)<.18)&(abs(yy)<.03))),
        'thin_spur_counterexample': (xx**2+yy**2 < .15**2)|((xx>.10)&(xx<.35)&(abs(yy)<.02)),
    }
    result = {name:raster_disk_coverage(a, (h,h)) for name,a in cases.items()}
    for name,a in cases.items():
        result[name]['topology_persistence'] = erosion_topology(a,(h,h))
    assert not result['thin_air_channel_width_0.06']['topology_persistence']['nominal_discrete_pass']
    assert not result['thin_neck_material_dumbbell']['topology_persistence']['nominal_discrete_pass']
    assert result['sharp_square_corners']['topology_persistence']['nominal_discrete_pass']
    assert result['thin_spur_counterexample']['topology_persistence']['nominal_discrete_pass']
    # These semantic checks must detect the intended thin phase regardless of
    # unrelated square corners elsewhere in the fixture.
    assert result['thin_island_width_0.06']['losses']['material']['removed_area'] > .01
    assert result['thin_air_channel_width_0.06']['losses']['air']['removed_area'] > .001
    assert result['thin_hole_diameter_0.06']['losses']['air']['removed_area'] > .002
    assert result['sharp_square_corners']['losses']['material']['removed_area'] > 0
    # Exact boundary tangency differs from a strict-clearance rule. Neither may
    # silently be substituted for the historical unspecified contact treatment.
    result['query_boundary_tangent_disk'] = dict(
        closed_containment=disk_in_rectangle((-.9,1.72), .12),
        strict_containment=disk_in_rectangle((-.9,1.72), .12, strict=True),
        note='Floating coordinates may resolve mathematical tangency within roundoff; use exact integer-scaled contract for enforcement.')
    points = np.array([[0.,1.6],[0.,1.6+1e-6]])
    result['query_boundary_material_contact'] = validate_material_samples(points,np.array([2.,1.]),np.array([0.,0.]))
    assert not result['query_boundary_material_contact']['query_air_pass']
    return result


def run_audit(old_root, output):
    started_cpu, started_wall = time.process_time(), time.monotonic()
    old_root, output = Path(old_root).resolve(), Path(output).resolve()
    decoder = types.ModuleType('frozen_decoder_only')
    decoder_path = old_root/'vendor/p0/geometry.py'
    # Execute this NumPy-only decoder without an import-loader bytecode write.
    exec(compile(decoder_path.read_text(), str(decoder_path), 'exec'), decoder.__dict__)
    config = json.loads((old_root/'docs/P0_FROZEN_CONFIG.json').read_text())
    if tuple(config['truth_seeds']) != SEEDS or config['materials']['minimum_feature'] != .1:
        raise ValueError('Frozen truth seeds/minimum_feature changed; audit contract requires review')
    records = {}
    h = .0025
    x, y = np.arange(-2+h/2, 2, h), np.arange(h/2, 3, h)
    xx, yy = np.meshgrid(x,y)
    points = np.c_[xx.ravel(), yy.ravel()]
    for seed in SEEDS:
        er, sigma = decoder.wall_values(points,seed)
        occupied = ((er != 1)|(sigma != 0)).reshape(xx.shape)
        assert np.array_equal(occupied.ravel(), wall_implicit(points,seed)>0)
        records[str(seed)] = dict(
            bounds_query_samples=validate_material_samples(points,er,sigma),
            analytic_envelope=dict(x=[-1.56,1.56],y=[.4446,1.4004],
                material_real_epsilon_on_occupied=[3.2,5.0], conductivity_on_occupied=[.018,.036],
                maximum_area_bound=3.12*.565, query_vertical_separation_lower_bound=.1996,
                status='ANALYTIC_PASS_BOUNDS_AREA_DOMAIN_QUERY_AIR'),
            area_ladder=[contour_area(seed,r) for r in (1,2,4)],
            cavity_boundary=boundary_audit(seed),
            rolling_disk_diagnostic=raster_disk_coverage(occupied,(h,h)),
            topology_persistence_ladder=[erosion_topology(occupied,(h,h))],
            face_connectivity_sensitivity=[erosion_topology(occupied,(h,h),connectivity=4)],
            original_minimum_feature_verdict='UNDETERMINED_UNSPECIFIED_CRITERION')
    # Fine full-domain topology includes the thin bridge in seed2; no crop that
    # could change phase connectivity. This is geometry only, no field arrays.
    hf=.00125
    xf,yf=np.arange(-2+hf/2,2,hf),np.arange(hf/2,3,hf)
    for seed in SEEDS:
        fine=np.empty((len(yf),len(xf)),dtype=bool)
        for j,yrow in enumerate(yf):
            fine[j]=wall_implicit(np.c_[xf,np.full(len(xf),yrow)],seed)>0
        records[str(seed)]['topology_persistence_ladder'].append(erosion_topology(fine,(hf,hf)))
        records[str(seed)]['face_connectivity_sensitivity'].append(erosion_topology(fine,(hf,hf),connectivity=4))
    probes=[]
    for center,radius,contrast in config['probes']:
        probes.append(dict(center=center,radius=radius,contrast=contrast,
                           query_containment=disk_in_rectangle(center,radius),
                           real_epsilon=1+contrast,conductivity=.02,
                           diameter=2*radius,added_area=float(np.pi*radius**2)))
    result=dict(status='ZERO_MAXWELL_AUDIT_COMPLETE_PRIOR_DEFINITION_HOLD',
        source_root=str(old_root),physical_solves_run=0,complex_rhs=0,lu=0,
        frozen_minimum_feature=config['materials']['minimum_feature'],
        minimum_feature_original_definition='scalar only: no width/ball/reach, phase, corner, boundary or tolerance contract',
        suggested_contract_status='NEW_ILLUSTRATIVE_ONLY_REQUIRES_EXPLICIT_ADOPTION',
        truth=records,probes=probes,
        all_probe_material_and_query_checks_pass=all(z['query_containment']['contained'] and 1 <= z['real_epsilon'] <=9 for z in probes),
        analytic_maximum_base_plus_one_probe_area=3.12*.565+np.pi*.12**2,
        adversarial_cases=adversarial_cases(),
        cpu_seconds=time.process_time()-started_cpu,wall_seconds=time.monotonic()-started_wall)
    output.mkdir(parents=True,exist_ok=True)
    (output/'COMPUTED.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old-root',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    r=run_audit(args.old_root,args.output)
    print(json.dumps({k:r[k] for k in ('status','physical_solves_run','complex_rhs','lu','cpu_seconds','wall_seconds')}))
