"""Controlled public fixed-support family, with an explicit freeze gate.

Import performs no file access, freeze, Book construction, data generation or
PDE solve. Only point(...) invokes a caller-funded Forward. No P1 truth path is
accepted or read. Nuisance gamma is external to the material/FEM coefficients.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import tempfile
import time
import types
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
FAMILY_ID='W_FIXED_SUPPORT_JOINT_V1'
PRIOR_PATH=Path('configs/W_CONTROLLED_PRIOR.json')
WALLS_PATH=Path('inputs/W_FAMILY_WALLS.npz')
SEED=2026102202
FREQUENCIES=(.9975,.99875,1.,1.00125,1.0025)
SIGMA_REAL=.01
U_BOUNDS=(-.01,.01)
CONTRAST_BOUNDS=(2.94,3.06)
CENTERS=((-0.9,2.15),(-0.8,2.15))
RADIUS=.12
TARGET_CONDUCTIVITY=.02


def _inside(root,path):
    p=Path(path)
    p=(root/p).resolve() if not p.is_absolute() else p.resolve()
    if not p.is_relative_to(root):raise ValueError('Path must remain within the controlled research directory')
    return p


def _read_json(path):
    return json.loads(Path(path).read_text())


def _write_exclusive(path,value):
    with Path(path).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.write('\n');stream.flush();os.fsync(stream.fileno())


def _q(value):
    q=np.asarray(value)
    if np.iscomplexobj(q) or q.shape!=(108,) or not np.isfinite(q).all() or np.any(abs(q)>1):
        raise ValueError('Each public wall requires108 finite real LinearWorld parameters in[-1,1]')
    return np.array(q,dtype=float,copy=True)


def _branch(value):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or int(value) not in (0,1):
        raise ValueError('branch must be integer0 or1')
    return int(value)


def _bounded(value,bounds,name):
    # Tolerance is exactly one float64 ULP outward from the SAME boundary.
    # It is not a wider physical interval or a relative/absolute tolerance.
    z=np.asarray(value)
    if z.ndim!=0 or np.iscomplexobj(z) or z.dtype.kind not in 'fiu':
        raise ValueError(name+' must be a finite real scalar')
    v=np.longdouble(z.item());lo,hi=(np.longdouble(t) for t in bounds)
    if not np.isfinite(v):raise ValueError(name+' must be finite')
    if v<lo:
        if v>=np.longdouble(np.nextafter(float(lo),-np.inf)):return float(lo)
        raise ValueError(name+' below frozen bound')
    if v>hi:
        if v<=np.longdouble(np.nextafter(float(hi),np.inf)):return float(hi)
        raise ValueError(name+' above frozen bound')
    return float(v)


def _acquisition(original):
    tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)]
    rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
    if (original['noise']['each_real_quadrature_std']!=SIGMA_REAL
        or tuple(original['frequencies']['narrow'])!=FREQUENCIES
        or not np.array_equal(np.asarray(original['tx']),tx)
        or not np.array_equal(np.asarray(original['rx']),rx)):
        raise RuntimeError('Original acquisition differs from the frozen controlled family')
    return dict(frequencies=list(FREQUENCIES),band='narrow',tx=tx.tolist(),rx=rx.tolist(),
        complex_response_shape=[5,4,12],sigma_each_real_quadrature=SIGMA_REAL,
        source_normalization=original['source_normalization'],
        whitening='full Re/Im norm divided by.01, no channel-count normalization')


def _validated_nomination(root,nomination):
    path=_inside(root,nomination);cfg=_read_json(path)
    if (cfg.get('seed')!=SEED or cfg.get('band')!='narrow' or cfg.get('kind')!='joint'
        or cfg.get('chart','linear')!='linear' or cfg.get('cavities',[])):
        raise RuntimeError('Family requires the fixed seed2/narrow/LinearWorld joint pair without cavity edits')
    identifier=cfg.get('id')
    if not isinstance(identifier,str) or Path(identifier).name!=identifier or identifier in ('.','..'):
        raise ValueError('Invalid nomination validation identifier')
    directory=root/'results/B2/validation'/identifier
    gate=_read_json(directory/'RESULTS.json')
    if gate.get('numerical_pass') is not True or gate.get('joint_physics_threshold') is not True:
        raise RuntimeError('Freeze blocked: independent numerical_pass AND joint_physics_threshold required')
    if _read_json(directory/'MANIFEST.json')!=cfg:
        raise RuntimeError('Nomination differs from independently validated manifest')
    distance=gate.get('d_joint_shifted')
    if (isinstance(distance,bool) or not isinstance(distance,(int,float)) or not np.isfinite(distance)
        or distance>1 or distance<0 or gate.get('target_shift')!=.1 or gate.get('mass_relative_change')!=0):
        raise RuntimeError('Validated joint distance/shift/mass fields are inconsistent with the frozen edit')
    wall_source=_inside(root,cfg['x_file'])
    with np.load(wall_source,allow_pickle=False) as archive:q1=_q(archive['x'])
    q0=np.zeros(108)
    for branch,q in ((0,q0),(1,q1)):
        endpoint=directory/f'fields/joint8_world{branch}_fine/endpoint_0000/PARAMETERS.json'
        recorded=_read_json(endpoint)
        if not np.array_equal(np.asarray(recorded['x']),q):
            raise RuntimeError('Current wall vector differs from independently validated fine endpoint')
        meta=recorded.get('metadata',{})
        if meta.get('seed')!=SEED or meta.get('cavities',[]):
            raise RuntimeError('Validated endpoint seed/support differs from controlled fixed-support pair')
        expected=(CENTERS[branch],RADIUS,3.)
        probe=meta.get('probe')
        if probe is None or not (np.array_equal(np.asarray(probe[0]),expected[0])
            and probe[1]==expected[1] and probe[2]==expected[2]):
            raise RuntimeError('Validated endpoint target differs from fixed branch center/radius/contrast')
    return cfg,q0,q1,path,gate,directory,wall_source


def freeze_family(nomination_path,*,root=None):
    """Gate, copy public q0/q1, then publish a NEW prior; no observations.

    Never overwrite either destination. A failed/incomplete publication is not
    automatically reused. This function must be invoked explicitly by the
    orchestrator after validation, never at import or Family.load time.
    """
    root=Path(root).resolve() if root is not None else ROOT
    prior,walls=root/PRIOR_PATH,root/WALLS_PATH
    if prior.exists() or walls.exists():raise FileExistsError('Controlled prior/walls already exist; refusing overwrite')
    cfg,q0,q1,nomination,gate,validation,source=_validated_nomination(root,nomination_path)
    original=_read_json(root/'configs/ORIGINAL_P0.json');acquisition=_acquisition(original)
    config=dict(family_id=FAMILY_ID,schema_version=1,frozen_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        identity='CONTROLLED_PUBLIC_FIXED_SUPPORT_TWO_WALL_SUBSET',
        intended_scope='independent restricted simulation family; not general unknown-wall online inference',
        original_minimum_feature_status='UNRESOLVED_ORIGINAL_SCALAR_NOT_REDEFINED',
        strict_original_prior_membership_claim=False,
        wall_support=dict(seed=SEED,definition='original analytic wall support shared by both branches',
            chart='LinearWorld54 nonnegative tensor hats/108 bounded affine parameters',
            wall_parameter_file=str(WALLS_PATH),npz_keys=['q0','q1'],wall_probability=[.5,.5],
            q0='zeros108',q1='fixed independently validated endpoint copied verbatim',cavities=[]),
        target=dict(branch_centers=[list(c) for c in CENTERS],radius=RADIUS,
            offset_x=dict(distribution='uniform',bounds=list(U_BOUNDS),independent_of_branch=True),
            contrast=dict(distribution='uniform',bounds=list(CONTRAST_BOUNDS),independent_of_branch_and_offset=True),
            conductivity=TARGET_CONDUCTIVITY,law='epsilon_target=1+contrast+i*.02/f',
            geometry='center=branch_center+(u,0); fixed radius.12'),
        nuisance=dict(shared='ONE gamma for every frequency/source/receiver; external multiply only',
            amplitude=dict(distribution='uniform',bounds=[.99,1.01]),
            phase_radians=dict(distribution='uniform',bounds=[-.01,.01]),independent_of_material=True),
        acquisition=acquisition,
        forward_numerics=dict(h=.13,tolerance=2.5e-5,contour_resolution=4,circle_segments=512,
            order=3,pml_strength=6,max_dofs=110000,max_factor_bytes=3500000000,stream_frequencies=True),
        provenance=dict(nomination=str(nomination.relative_to(root)),nomination_snapshot=cfg,
            independent_validation=str(validation.relative_to(root)),validation_gate_snapshot=gate,
            wall_parameter_source=str(source.relative_to(root)),
            publication_gate='independent numerical_pass AND joint_physics_threshold',
            conditional_fixed_support_scope=True),
        observation_or_truth_generation='NONE; handled only by a separate explicitly authorized process',
        probability_origin='new explicit controlled prior masses; not inferred from witness existence or likelihood fit',
        boundary_roundoff='only one outward float64ULP maps back to its same frozen boundary')
    prior.parent.mkdir(parents=True,exist_ok=True);walls.parent.mkdir(parents=True,exist_ok=True)
    # Temporary files plus hard-link publication give atomic NO-REPLACE writes.
    # The JSON is published last and acts as the complete-family marker.
    with tempfile.TemporaryDirectory(prefix='.family-freeze-',dir=root) as temporary:
        temp=Path(temporary);wall_temp=temp/'walls.npz';prior_temp=temp/'prior.json'
        with wall_temp.open('xb') as stream:
            np.savez_compressed(stream,q0=q0,q1=q1)
            stream.flush();os.fsync(stream.fileno())
        _write_exclusive(prior_temp,config)
        if prior.exists() or walls.exists():raise FileExistsError('Concurrent existing family detected; refusing overwrite')
        os.link(wall_temp,walls)
        try:os.link(prior_temp,prior)
        except BaseException:
            # Own new wall link only; no prior was published by this operation.
            # Keep it as explicit incomplete-publication evidence, never overwrite.
            raise RuntimeError('Prior publication failed after wall copy; incomplete wall file retained for review')
    return dict(family_id=FAMILY_ID,prior=str(prior),walls=str(walls),
        observations_generated=0,truth_generated=0,physical_solves_run=0)


class Family:
    """Public controlled prior loader and explicit, caller-funded forward map."""
    def __init__(self,root,config,walls):
        self.root=Path(root).resolve();self.config=config
        self.walls=tuple(_q(q) for q in walls)
        for q in self.walls:q.flags.writeable=False

    @classmethod
    def load(cls,root=None):
        root=Path(root).resolve() if root is not None else ROOT
        config=_read_json(root/PRIOR_PATH)
        if config.get('family_id')!=FAMILY_ID:raise RuntimeError('Unexpected public controlled family')
        wall=config.get('wall_support',{});target=config.get('target',{});nuisance=config.get('nuisance',{})
        if (wall.get('seed')!=SEED or wall.get('wall_parameter_file')!=str(WALLS_PATH)
            or wall.get('wall_probability')!=[.5,.5] or wall.get('cavities')!=[]
            or target.get('branch_centers')!=[list(c) for c in CENTERS]
            or target.get('radius')!=RADIUS or target.get('conductivity')!=TARGET_CONDUCTIVITY
            or target.get('offset_x',{}).get('bounds')!=list(U_BOUNDS)
            or target.get('contrast',{}).get('bounds')!=list(CONTRAST_BOUNDS)
            or nuisance.get('amplitude',{}).get('bounds')!=[.99,1.01]
            or nuisance.get('phase_radians',{}).get('bounds')!=[-.01,.01]
            or any(z.get('distribution')!='uniform' for z in (target.get('offset_x',{}),target.get('contrast',{}),nuisance.get('amplitude',{}),nuisance.get('phase_radians',{})))):
            raise RuntimeError('Public prior differs from fixed controlled family contract')
        expected_numerics=dict(h=.13,tolerance=2.5e-5,contour_resolution=4,circle_segments=512,
            order=3,pml_strength=6,max_dofs=110000,max_factor_bytes=3500000000,stream_frequencies=True)
        if config.get('forward_numerics')!=expected_numerics:
            raise RuntimeError('Public prior forward numerics differ from fixed evaluator')
        acquired=config.get('acquisition',{})
        original=dict(noise={'each_real_quadrature_std':acquired.get('sigma_each_real_quadrature')},
            frequencies={'narrow':acquired.get('frequencies',[])},tx=acquired.get('tx',[]),rx=acquired.get('rx',[]),
            source_normalization=acquired.get('source_normalization'))
        _acquisition(original)
        gates=config.get('provenance',{}).get('validation_gate_snapshot',{})
        if gates.get('numerical_pass') is not True or gates.get('joint_physics_threshold') is not True:
            raise RuntimeError('Frozen prior lacks required independent validation gate snapshot')
        with np.load(root/WALLS_PATH,allow_pickle=False) as archive:
            walls=(_q(archive['q0']),_q(archive['q1']))
        if not np.array_equal(walls[0],np.zeros(108)):raise RuntimeError('Public branch0 must be zeros108')
        return cls(root,config,walls)

    def coordinates(self,b,u,c):
        b=_branch(b);u=_bounded(u,U_BOUNDS,'u');c=_bounded(c,CONTRAST_BOUNDS,'contrast')
        center=(CENTERS[b][0]+u,CENTERS[b][1])
        return dict(branch=b,u=u,contrast=c,center=center,radius=RADIUS,
            conductivity=TARGET_CONDUCTIVITY,branch_prior_mass=.5,
            coordinate_role='PUBLIC_PRIOR_NODE_NOT_GT_LABEL',nuisance_in_material=False)

    def material(self,b,u,c,grid,f=1.):
        """Decode the public material, with no FEM, Book or nuisance multiplier."""
        node=self.coordinates(b,u,c)
        if not np.isscalar(f) or np.iscomplexobj(f) or not np.isfinite(f) or f<=0:raise ValueError('frequency must be finite positive real')
        z=np.asarray(grid)
        if np.iscomplexobj(z) or z.ndim<2 or z.shape[-1]!=2 or not np.isfinite(z).all():
            raise ValueError('grid must have finite real shape(...,2)')
        from p0.geometry import wall_values
        from .endpoint_geometry import linear_material_values
        decoder=types.SimpleNamespace(wall_values=wall_values)
        probe=(node['center'],RADIUS,node['contrast'])
        er,sigma=linear_material_values(SEED,self.walls[node['branch']],z.reshape(-1,2),decoder,probe=probe)
        return (er+1j*np.divide(sigma,float(f))).reshape(z.shape[:-1])

    def point(self,b,u,c,book,folder):
        """Run one unique complete complex F endpoint only when explicitly called.

        Book is provided by the orchestrator; no book/budget is created here.
        Every frequency is immediately saved by Forward. No Jacobian, no gain
        fitting, no noise/truth generation, no cached LU across endpoints.
        """
        node=self.coordinates(b,u,c)
        if self.root!=ROOT:
            raise RuntimeError('Existing Forward is bound to this campaign ROOT; other roots are load/material-only')
        if not all(callable(getattr(book,k,None)) for k in ('reserve','event','totals')):
            raise TypeError('A caller-authorized Book-like reservation ledger is required')
        if hasattr(book,'root') and Path(book.root).resolve()!=self.root:
            raise ValueError('Book belongs to another research directory')
        results=self.root/'results';path=Path(folder)
        dest=_inside(self.root,results/path if not path.is_absolute() else path)
        if not dest.is_relative_to(results):raise ValueError('point folder must be inside this family results directory')
        dest.mkdir(parents=True,exist_ok=False)
        _write_exclusive(dest/'PUBLIC_PRIOR_NODE.json',dict(family_id=FAMILY_ID,coordinates=node,
            public_prior=str(PRIOR_PATH),wall_parameters=str(WALLS_PATH),status='PREPARED_NOT_SOLVED',
            full_response_shape=[5,4,12],nuisance_application='external multiply only'))
        # Physical dependencies are intentionally lazy: importing/load/material
        # cannot construct a solver. Forward itself enforces reservations.
        from .native_world import LinearWorld
        from .physical_witness import mesh_for
        from .physics import Forward
        outer=self
        class PublicPriorWorld(LinearWorld):
            def material_metadata(model,x):
                result=super().material_metadata(x)
                result.update(family_id=FAMILY_ID,identity='PUBLIC_CONTROLLED_PRIOR_NODE',
                    prior_coordinates=node,public_prior=str(PRIOR_PATH),
                    public_wall_parameter_file=str(WALLS_PATH),
                    conditional_fixed_support_family=True,original_minimum_feature_status='UNRESOLVED',
                    coordinate_role='PUBLIC_PRIOR_NODE_NOT_GT_LABEL',nuisance_in_coefficient=False)
                return result
        try:
            probe=(node['center'],RADIUS,node['contrast'])
            model=PublicPriorWorld(SEED,probe=probe,cavities=())
            mesh=mesh_for(SEED,probe,h=.13,tol=2.5e-5,res=4,n=512,cavities=())
            forward=Forward(model,mesh,FREQUENCIES,book,str(dest.relative_to(results)),
                phase='P1_controlled_public_family_forward',identity='PUBLIC_CONTROLLED_PRIOR_NODE',
                h=.13,order=3,pml_strength=6,max_dofs=110000,max_factor_bytes=3500000000,
                retain_current=False)
            data,_=forward.evaluate(outer.walls[node['branch']].copy(),jac=False)
            data=np.asarray(data)
            if data.shape!=(5,4,12) or not np.iscomplexobj(data) or not np.isfinite(data).all():
                raise RuntimeError('Full raw complex response must have finite shape(5,4,12)')
            with (dest/'FULL_RESPONSE.npz').open('xb') as stream:
                np.savez_compressed(stream,F=data,branch=node['branch'],u=node['u'],contrast=node['contrast'],
                    frequencies=np.asarray(FREQUENCIES))
            _write_exclusive(dest/'POINT_COMPLETE.json',dict(family_id=FAMILY_ID,coordinates=node,
                complete_complex_response_shape=[5,4,12],frequency_count=5,
                nuisance_applied=False,noise_generated=False,truth_generated=False))
            return data.copy()
        except BaseException as exc:
            _write_exclusive(dest/'POINT_FAILURE.json',dict(error_type=type(exc).__name__,error=str(exc),
                coordinates=node,completed_frequency_files='see endpoint_0000; preserved without recovery solve'))
            raise
