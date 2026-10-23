"""Fixed-witness parameter/position sensitivity and compensation ablations.

No optimization or added gain. Perturbations are stress diagnostics, not
additional qualified worlds (a negative position perturbation is below the
original .1 target-shift gate). All new actions use B2 validation allowance.
"""
import sys,json
import numpy as np
from .accounting import ROOT,ScopedBook,write
from .physical_witness import mesh_for
from .native_world import LinearWorld
from .physics import Forward

def main(manifest):
    cfg=json.loads((ROOT/manifest).read_text());vdir=ROOT/'results/B2/validation'/cfg['id']
    gate=json.loads((vdir/'RESULTS.json').read_text())
    if not gate['joint_physics_threshold'] or not gate['numerical_pass']:raise RuntimeError('no qualified joint endpoint')
    book=ScopedBook('B2','validation');out=vdir/'stress';out.mkdir(exist_ok=False)
    freq=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())['frequencies'][cfg['band']]
    q=np.load(ROOT/cfg['x_file'])['x'];reference=np.load(vdir/'joint8_PAIR.npz')['fine0']
    original=np.load(vdir/'joint8_PAIR.npz')['fine1'];seed=cfg['seed'];rows=[]
    cases=[('parameter_plus',np.clip(q+.002,-1,1),0.),('parameter_minus',np.clip(q-.002,-1,1),0.),
        ('position_plus',q,.005),('position_minus',q,-.005),
        ('real_material_only',np.r_[q[:54],np.zeros(54)],0.),
        ('conductivity_only',np.r_[np.zeros(54),q[54:]],0.)]
    for name,point,dx in cases:
        probe=((-0.8+dx,2.15),.12,3.);model=LinearWorld(seed,probe)
        fw=Forward(model,mesh_for(seed,probe,h=.13,tol=2.5e-5,res=4,n=512),freq,book,
            str((out/name).relative_to(ROOT/'results')),phase='B2_fixed_witness_stress',identity='ORACLE_STRESS_NOT_NEW_DISCOVERY',
            h=.13,pml_strength=6,max_dofs=110000,max_factor_bytes=3500000000)
        d,_=fw.evaluate(point)
        row=dict(case=name,d_joint=float(np.linalg.norm(d-reference)/.01),
            change_from_nominal_world2=float(np.linalg.norm(d-original)/.01),
            target_shift=.1+dx,remaining_under_d1=bool(np.linalg.norm(d-reference)/.01<=1),
            new_qualified_witness_claim=False,compensation_not_reoptimized=True)
        rows.append(row);write(out/'RESULTS.json',rows)
        np.savez_compressed(out/f'{name}_PAIR.npz',q=point,data=d,reference=reference,nominal_world2=original)
    write(out/'COMPLETE.json',dict(totals=book.totals(),new_rhs=120))

if __name__=='__main__':main(sys.argv[1])
