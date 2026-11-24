"""Prospective spatial branch menu scored on two exposed native half-ring observations."""

import argparse
import json

import numpy as np
from flint import arb
from threadpoolctl import threadpool_limits

from branches import propose_branches
from native_box import HERE, OPTICS, NativeBox, ball, clock, elapsed, frobenius_squared, hi
from restricted import score_hard

SLIT='optfresh_01__pw180_10x20__n02'
SHELL='s2d_004__pw180_10x20'


def inputs(name):
    if name==SLIT:
        source=OPTICS/f'synthetic/public/{name}.npz'
        with np.load(source,allow_pickle=False) as f:
            y=f['Y']
        audit=json.loads((OPTICS/'results/theory/fine_feasibility.json').read_text())['observations']
        record=next(r for r in audit if r['id']==name)
        return y,record['inflated_radius'],str(source), 'verified native generating material lies in this padded ball'
    source=OPTICS.parent/f'srap_cie_array_study_20261121/data/public/{name}.npz'
    with np.load(source,allow_pickle=False) as f:
        y,delta=f['Y'],float(f['delta'])
        meta=json.loads(str(f['case_json']))
    assert (meta['fine'],meta['frequency'],meta['array'],meta['noise_relative'])==(96,400e6,'PW180_10x20',.1)
    radius=hi(arb(delta)+frobenius_squared(ball(y)).sqrt()/10**6)
    return y,radius,str(source),'declared sample radius plus 1e-6 relative frozen-forward allowance'


def run(args):
    start=clock()
    hw_start=clock()
    with np.load(OPTICS/'synthetic/public/fine_geometry.npz',allow_pickle=False) as f:
        cert=NativeBox(f['spectrum'],f['PW180_10x20__S'],f['PW180_10x20__B'])
    prepared=elapsed(hw_start)
    records={}
    for name in args.case or [SLIT,SHELL]:
        initial=clock()
        y,radius,source,radius_scope=inputs(name)
        fold=HERE.parent/'diagnosis/trajectories'/name
        saved=fold/'step_0064.npz'
        with np.load(saved,allow_pickle=False) as f:
            proposed=propose_branches(f['m'],f['J'],tile_count=args.tiles,foreground=args.upper,background=args.background)
        with np.load(fold/'materials.npz',allow_pickle=False) as f:
            idx=np.flatnonzero(f['step']==64)
            prefix_wall=float(f['engine_wall'][idx[0]])
        preparation=elapsed(initial)
        row=dict(observation=name,source=str(source),radius=radius,radius_scope=radius_scope,
                 observation_norm=float(np.linalg.norm(y)),prefix_source=str(saved.relative_to(HERE.parent)),
                 prefix_material_side=32,prefix_solver_wall_seconds=prefix_wall,preparation=preparation,
                 tile_score=proposed['tile_score'],branches=[])
        records[name]=row
        args.output.write_text(json.dumps(dict(preparation=prepared,by_case=records,**elapsed(start)),indent=2)+'\n')
        for branch in proposed['branches']+proposed['development_hard_controls']:
            if args.branch and branch.name not in args.branch:
                continue
            if branch.name=='wildcard':
                result=dict(branch=branch.name,decision='RETAIN',rigor='always_present_open_fallback',radius=radius,
                            residual_lower=0.,wall_seconds=0.,cpu_seconds=0.)
            elif branch.name=='uniform_low_050':
                result=cert.score(y,radius,branch,t_grid=(.8,0.,1.),c_grid=(.1875,.125,.0625),
                                  alpha_grid=(.0005,.00025,.001,2.**-13),max_proofs=1)
            elif args.active and 'hard_support' in branch.name:
                result=score_hard(cert,y,radius,branch,max_proofs=1)
            else:
                result=cert.score(y,radius,branch,t_grid=(0.,.5,1.), c_grid=(.9375,.5,.375),max_proofs=1)
            result['tile_caps']=branch.caps
            row['branches'].append(result)
            args.output.write_text(json.dumps(dict(preparation=prepared,by_case=records,**elapsed(start)),indent=2)+'\n')
            best=max((p.get('numeric_residual',0.) for p in result.get('screen',[])),default=0.)
            nus=[p['nu'] for p in result.get('screen',[]) if 'nu' in p]
            print(name,branch.name,'numeric',best,'minimum_nu',min(nus) if nus else None,
                  'strict',result['residual_lower'],result['decision'],result['wall_seconds'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',action='append',choices=(SLIT,SHELL))
    parser.add_argument('--branch',action='append')
    parser.add_argument('--tiles',type=int,default=4)
    parser.add_argument('--upper',type=float,default=4.)
    parser.add_argument('--background',type=float,default=.05)
    parser.add_argument('--active',action='store_true')
    parser.add_argument('--output',type=type(HERE),default=HERE/'development.json')
    with threadpool_limits(1):
        run(parser.parse_args())
