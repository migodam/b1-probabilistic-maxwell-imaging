"""Pure decoder/geometry checks for LinearWorld endpoints; no physics imports.

This module does not import native_world: its parent imports physics machinery.
It implements the published 6x9 partition-of-unity chart directly for checks.
Original minimum_feature=.1 semantics remain unresolved. Topology diagnostics
are explicitly weaker than a lower bound on every local width.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time
import types
import numpy as np
from scipy.optimize import differential_evolution
from .geometry_audit import (SEEDS, QUERY, DOMAIN, phases, wall_implicit,
    cavity_boundary, strip_boundary, boundary_audit, erosion_topology,
    disk_in_rectangle, validate_material_samples)

FIXED_CAVITY = ((-.85,.85),.06)


def trig_range(a,b,*,cosine=False):
    if cosine: a,b=a+np.pi/2,b+np.pi/2
    if b-a >= 2*np.pi: return -1.,1.
    v=[np.sin(a),np.sin(b)]
    for k in range(int(np.ceil((a-np.pi/2)/np.pi)),int(np.floor((b-np.pi/2)/np.pi))+1):
        v.append(np.sin(np.pi/2+k*np.pi))
    return float(min(v)),float(max(v))


def cavity_clearance(seed,center=FIXED_CAVITY[0],radius=FIXED_CAVITY[1]):
    """Analytic interval containment proof plus numerical Euclidean distances."""
    c=np.asarray(center,dtype=float);r=float(radius)
    if c.shape!=(2,) or not np.isfinite(c).all() or not np.isfinite(r) or r<=0:
        raise ValueError('Finite center(2), positive radius required')
    p=phases(seed);xl,xr=c[0]-r,c[0]+r;yb,yt=c[1]-r,c[1]+r
    smin,smax=trig_range(2.4*xl+p[0],2.4*xr+p[0])
    low_min=.64+.09*xl+.055*smin;low_max=.64+.09*xr+.055*smax
    tmin,tmax=trig_range(2.1*xl+p[1],2.1*xr+p[1],cosine=True)
    top_min=low_min+.48+.085*tmin
    smn,_=trig_range(4*yb+p[2],4*yt+p[2])
    side_margin=1.48+.08*smn-max(abs(xl),abs(xr))
    # Reverse triangle inequality in the ellipse-scaled norm bounds the entire
    # disk away from the old cavity, whose radial right side never exceeds1.18.
    scaled_center=float(np.hypot((c[0]+.2)/.25,(c[1]-.92)/.15))
    cavity_implicit_lower=scaled_center-r/.15-1.18
    interval=dict(bottom_vertical_margin_lower=float(yb-low_max),
        top_vertical_margin_lower=float(top_min-yt),side_horizontal_margin_lower=float(side_margin),
        old_cavity_implicit_margin_lower=float(cavity_implicit_lower),
        method='interval trig extrema and reverse triangle inequality; conservative analytic bounds')
    inside=all(interval[k]>0 for k in ('bottom_vertical_margin_lower','top_vertical_margin_lower',
                                      'side_horizontal_margin_lower','old_cavity_implicit_margin_lower'))
    distances={}
    for side in ('bottom','top'):
        opt=differential_evolution(lambda v:float((v[0]-c[0])**2+
            (strip_boundary(v[0],seed,side)-c[1])**2),[(-2.,2.)],seed=0,tol=1e-12,polish=True)
        distances[side]=dict(center_to_curve=float(np.sqrt(max(0,opt.fun))),
            disk_boundary_clearance=float(np.sqrt(max(0,opt.fun))-r),
            nearest_point=[float(opt.x[0]),float(strip_boundary(opt.x[0],seed,side))])
    for sign,name in ((-1,'left'),(1,'right')):
        def side_point(y): return np.array([sign*(1.48+.08*np.sin(4*y+p[2])),y])
        opt=differential_evolution(lambda v:float(np.sum((side_point(v[0])-c)**2)),
            [(0.,3.)],seed=0,tol=1e-12,polish=True)
        distances[name]=dict(center_to_curve=float(np.sqrt(max(0,opt.fun))),
            disk_boundary_clearance=float(np.sqrt(max(0,opt.fun))-r),nearest_point=side_point(opt.x[0]).tolist())
    opt=differential_evolution(lambda v:float(np.sum((cavity_boundary(v[0],seed)-c)**2)),
        [(0.,2*np.pi)],seed=0,tol=1e-12,polish=True)
    distances['original_cavity']=dict(center_to_curve=float(np.sqrt(max(0,opt.fun))),
        disk_boundary_clearance=float(np.sqrt(max(0,opt.fun))-r),nearest_point=cavity_boundary(opt.x[0],seed).tolist())
    return dict(seed=int(seed),center=c.tolist(),radius=r,diameter=2*r,
        entire_closed_disk_in_original_occupied_wall_analytic=bool(inside),
        analytic_interval_bounds=interval,euclidean_clearances=distances,
        euclidean_distance_status='numerical minimization; no certified global rounding error bound',
        removed_area_exact_if_contained=float(np.pi*r*r) if inside else None,
        feature_width_warning='Diameter alone does not establish material bridge width. Original feature semantics unresolved.')


def _chart_values(q,points):
    """Direct bilinear interpolation of the 6x9 nonnegative hat coefficients."""
    z=np.asarray(points,dtype=float)
    x=np.clip((z[:,0]+1.6)/.4,0.,8.);y=np.clip((z[:,1]-.4)/.2,0.,5.)
    ix=np.minimum(x.astype(int),7);iy=np.minimum(y.astype(int),4)
    ax=x-ix;ay=y-iy;a=q.reshape(6,9)
    return ((1-ay)*((1-ax)*a[iy,ix]+ax*a[iy,ix+1])
            +ay*((1-ax)*a[iy+1,ix]+ax*a[iy+1,ix+1]))


def linear_material_values(seed,parameters,points,decoder,*,cavities=(),probe=None):
    """Pointwise independent decoder; does not import the physical chart module."""
    q=np.asarray(parameters)
    if np.iscomplexobj(q) or q.shape!=(108,) or not np.isfinite(q).all() or np.any(abs(q)>1):
        raise ValueError('LinearWorld requires108 finite real parameters in[-1,1]')
    z=np.asarray(points,dtype=float)
    if z.ndim!=2 or z.shape[1]!=2 or not np.isfinite(z).all():
        raise ValueError('points must be finite(N,2)')
    er,sigma=decoder.wall_values(z,seed)
    for center,radius in cavities:
        active=np.linalg.norm(z-np.asarray(center),axis=1)<radius
        er[active]=1.;sigma[active]=0.
    er=1+(er-1)*(1+.6*_chart_values(q[:54],z))
    sigma=sigma*(1+.75*_chart_values(q[54:],z))
    if probe is not None:
        center,radius,contrast=probe
        active=np.linalg.norm(z-np.asarray(center),axis=1)<radius
        er[active]=1+contrast;sigma[active]=.02
    return er,sigma


def check_linear_world_endpoint(seed,parameters,*,cavities=(),probe=None):
    """Executable zero-physics endpoint contract and explicit unresolved flags.

    Analytic chart bounds are global, not inferred from a raster. Caveats on
    original features and seed1 notch are preserved even for valid q. The
    known-reference support chart is an ORACLE chart, not online inference.
    Geometry diagnostics may be computed separately at chosen resolutions.
    """
    q=np.asarray(parameters)
    if np.iscomplexobj(q) or q.shape!=(108,) or not np.isfinite(q).all() or np.any(abs(q)>1):
        raise ValueError('LinearWorld requires108 finite real parameters in[-1,1]')
    if seed not in SEEDS: raise ValueError('Only original frozen seeds audited')
    amin,amax=float(q[:54].min()),float(q[:54].max())
    smin,smax=float(q[54:].min()),float(q[54:].max())
    contrast_bounds=[2.2*(1+.6*amin),4.*(1+.6*amax)]
    conductivity_bounds=[.018*(1+.75*smin),.036*(1+.75*smax)]
    cavity_checks=[cavity_clearance(seed,c,r) for c,r in cavities]
    pair_disjoint=all(np.linalg.norm(np.asarray(cavities[i][0])-np.asarray(cavities[j][0]))>
        cavities[i][1]+cavities[j][1] for i in range(len(cavities)) for j in range(i))
    subtract=sum(np.pi*r*r for c,r in cavities) if pair_disjoint and all(
        z['entire_closed_disk_in_original_occupied_wall_analytic'] for z in cavity_checks) else 0.
    area_bound=1.7628-subtract
    probe_check=None
    if probe is not None:
        center,radius,contrast=probe
        if not np.isfinite(contrast) or not 0<=contrast<=8:
            raise ValueError('probe contrast must lie in[0,8]')
        probe_check=disk_in_rectangle(center,radius,QUERY)
        if not probe_check['contained']:
            raise ValueError('Probe must be contained in original queryT')
        area_bound+=np.pi*radius**2
    global_contrast_max=max(contrast_bounds[1],probe[2] if probe is not None else 0.)
    global_sigma_max=max(conductivity_bounds[1],.02 if probe is not None else 0.)
    original_boundary=boundary_audit(seed)
    return dict(model='LinearWorld6x9_hat_affine',seed=int(seed),parameter_box_pass=True,
        identity='ORACLE_KNOWN_REFERENCE_CHART_NOT_ONLINE',
        global_material_bounds_pass=bool(global_contrast_max<=8 and global_sigma_max<=.15),
        global_passivity_pass=bool(conductivity_bounds[0]>=0),
        positive_modulation_preserves_fixed_support=True,
        occupied_wall_contrast_interval=contrast_bounds,
        occupied_wall_real_epsilon_interval=[1+contrast_bounds[0],1+contrast_bounds[1]],
        occupied_wall_conductivity_interval=conductivity_bounds,
        all_domain_real_epsilon_interval=[1.,1+global_contrast_max],
        all_domain_conductivity_interval=[0.,global_sigma_max],
        bounds_proof='nonnegative tensor hats sum1 globally after clamping; q in[-1,1]',
        common_cross_frequency_law='epsilon_infinity+i conductivity/f, same two real fields',
        original_support_envelope=dict(x=[-1.56,1.56],y=[.4446,1.4004]),
        calibration_query_air_pass=probe is None,
        task_probe_query_containment=probe_check,
        query_air_note='A task probe intentionally occupiesT; calibration air condition applies to no-probe state.',
        maximum_occupied_area_analytic_upper=float(area_bound),area_limit_pass=bool(area_bound<=3),
        cavity_checks=cavity_checks,cavity_pairwise_disjoint=bool(pair_disjoint),
        original_seed_geometry=original_boundary,
        seed1_original_notch_preserved=bool(seed==SEEDS[0] and all(
            z['entire_closed_disk_in_original_occupied_wall_analytic'] for z in cavity_checks)),
        original_minimum_feature='UNRESOLVED_SCALAR_SEMANTICS',
        local_width_certified=False,physical_solves_run=0,complex_rhs=0,lu=0)


def run(old_root,root):
    cpu,wall=time.process_time(),time.monotonic()
    old_root,root=Path(old_root).resolve(),Path(root).resolve()
    decoder=types.ModuleType('frozen_decoder_no_bytecode')
    path=old_root/'vendor/p0/geometry.py'
    exec(compile(path.read_text(),str(path),'exec'),decoder.__dict__)
    baseline=json.loads((root/'AUDIT/geometry/COMPUTED.json').read_text())
    probe=decoder.probes()[8]
    shifted=((probe[0][0]+.1,probe[0][1]),probe[1],probe[2])
    probes=dict(original=probe,shifted=shifted,query_original=disk_in_rectangle(probe[0],probe[1]),
        query_shifted=disk_in_rectangle(shifted[0],shifted[1]),
        centroid_translation_exact=[.1,0.],dielectric_mass_relative_change_exact=0.,
        dielectric_mass_exact=float(probe[2]*np.pi*probe[1]**2),
        wall_separation_vertical_lower=float(shifted[0][1]-shifted[1]-1.4004))
    # Pointwise checker verification, without importing the solver-bearing
    # LinearWorld parent module. Cavity/probe overrides must be exact.
    fixtures=np.array([[-.85,.85],[-.9,.75],[-.2,.92],shifted[0]])
    checks={}
    for value in (-1.,0.,1.):
        er,si=linear_material_values(SEEDS[1],np.full(108,value),fixtures,decoder,
            cavities=(FIXED_CAVITY,),probe=shifted)
        assert er[0]==1. and si[0]==0.
        assert er[-1]==1+shifted[2] and si[-1]==.02
        assert np.all((er>=1)&(er<=9)&(si>=0)&(si<=.15))
        checks[str(value)]=dict(real_epsilon=er.tolist(),conductivity=si.tolist())
    try:
        check_linear_world_endpoint(SEEDS[1],np.full(108,1.001))
    except ValueError:
        checks['out_of_box_rejected']=True
    else:
        raise AssertionError('out-of-box vector accepted')
    records={}
    q=np.zeros(108)
    for seed in SEEDS:
        ep=check_linear_world_endpoint(seed,q,cavities=(FIXED_CAVITY,))
        topologies=[]
        for h in (.0025,.00125):
            x,y=np.arange(-2+h/2,2,h),np.arange(h/2,3,h)
            mask=np.empty((len(y),len(x)),dtype=bool)
            for j,yy in enumerate(y):
                pts=np.c_[x,np.full(len(x),yy)]
                er,si=linear_material_values(seed,q,pts,decoder,cavities=(FIXED_CAVITY,))
                mask[j]=(er!=1)|(si!=0)
            topologies.append(dict(h=h,connectivity8=erosion_topology(mask,(h,h)),
                connectivity4=erosion_topology(mask,(h,h),connectivity=4),sampled_area=float(mask.sum()*h*h)))
        areas=baseline['truth'][str(seed)]['area_ladder']
        ep['base_area_approximations']=areas
        ep['new_area_approximations']=[dict(resolution=a['resolution'],area=a['occupied_area']-np.pi*.06**2,
            status='old contour approximation minus analytically contained circle area') for a in areas]
        ep['new_cavity_topology_diagnostics']=topologies
        ep['shifted_probe_endpoint']=check_linear_world_endpoint(seed,q,probe=shifted)
        ep['parameter_box_extrema_proof']={name:check_linear_world_endpoint(seed,np.full(108,value))[
            'all_domain_real_epsilon_interval'] for name,value in [('all_minus1',-1.),('all_plus1',1.)]}
        records[str(seed)]=ep
    result=dict(status='ZERO_PHYSICS_FIXED_B2_EDIT_GEOMETRY_AUDIT',fixed_cavity=FIXED_CAVITY,
        probes=probes,seed_endpoints=records,decoder_self_checks=checks,physical_solves_run=0,complex_rhs=0,lu=0,
        original_minimum_feature='UNRESOLVED_SCALAR_SEMANTICS',
        warning='Topology persistence is not all-local-width certification; seed1 new top bridge clearance<.1.',
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.monotonic()-wall)
    out=root/'AUDIT/geometry/B2_EDITS.json'
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old-root',required=True);p.add_argument('--root',required=True)
    a=p.parse_args();r=run(a.old_root,a.root)
    print(json.dumps({k:r[k] for k in ('status','physical_solves_run','complex_rhs','lu','cpu_seconds','wall_seconds')}))
