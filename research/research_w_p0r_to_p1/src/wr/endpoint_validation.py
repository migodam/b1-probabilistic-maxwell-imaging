"""Independent finer, interface-fitted evaluation of a NOMINATED fixed pair.

No reoptimization occurs in this validator. Its result is an empirical
two-resolution check, not a rigorous continuous-PDE error enclosure.
"""
import sys,json
import numpy as np
from p0.geometry import probes
from .native_world import LinearWorld
from .physical_witness import ModulatedWorld,mesh_for
from .physics import Forward
from .accounting import ROOT,ScopedBook,write

def main(manifest):
    cfg=json.loads((ROOT/manifest).read_text());seed=cfg['seed'];band=cfg['band'];kind=cfg['kind']
    book=ScopedBook('B2','validation');out=ROOT/'results/B2/validation'/cfg['id']
    out.mkdir(parents=True,exist_ok=False)
    freq=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())['frequencies'][band]
    cls=LinearWorld if cfg.get('chart','linear')=='linear' else ModulatedWorld
    q=np.load(ROOT/cfg['x_file'])['x'];zero=np.zeros_like(q)
    cavities=tuple((tuple(x[0]),x[1]) for x in cfg.get('cavities',[]))
    original=probes()[8]
    shifted=((original[0][0]+.1,original[0][1]),original[1],original[2])
    write(out/'MANIFEST.json',cfg)
    def evaluate(which,probe=None,fine=True,label=''):
        model=cls(seed,probe=probe,cavities=cavities if which else ())
        h=.13 if fine else .145;tol=2.5e-5 if fine else 5e-5;res=4 if fine else 2;n=512 if fine else 256
        fw=Forward(model,mesh_for(seed,probe,h=h,tol=tol,res=res,n=n,cavities=cavities if which else ()),
            freq,book,str((out/'fields'/label).relative_to(ROOT/'results')),
            phase='B2_independent_endpoint_validation',identity='ORACLE_PHYSICAL_WITNESS_VALIDATION',
            h=h,pml_strength=6,max_dofs=110000,max_factor_bytes=3500000000)
        return fw.evaluate(q if which else zero)[0]
    rows=[]
    # Validate calibration and the directly claimed joint/probe8 endpoint on
    # both grids, including the correct circular interface for every target.
    for mode in ['cal','joint8']:
        pr0=None if mode=='cal' else original
        pr1=None if mode=='cal' else (shifted if kind=='joint' else original)
        co=[];fi=[]
        for which,pr in [(0,pr0),(1,pr1)]:
            co.append(evaluate(which,pr,False,f'{mode}_world{which}_coarse'))
            fi.append(evaluate(which,pr,True,f'{mode}_world{which}_fine'))
        row=dict(mode=mode,coarse_distance=float(np.linalg.norm(co[1]-co[0])/.01),
                 fine_distance=float(np.linalg.norm(fi[1]-fi[0])/.01),
                 world0_grid_distance=float(np.linalg.norm(co[0]-fi[0])/.01),
                 world1_grid_distance=float(np.linalg.norm(co[1]-fi[1])/.01))
        rows.append(row);write(out/'GRID_COMPARISON.json',rows)
        np.savez_compressed(out/f'{mode}_PAIR.npz',coarse0=co[0],coarse1=co[1],fine0=fi[0],fine1=fi[1])
    # Original d_T is a supremum over the fixed task class, not the norm of
    # 24 concatenated experiments. Full responses include wall scattering.
    tasks=[];arrays=[]
    for index,probe in enumerate(probes()):
        if index==8 and kind!='joint':
            pair=np.load(out/'joint8_PAIR.npz');d0,d1=pair['fine0'],pair['fine1']
        else:
            d0=evaluate(0,probe,True,f'probe{index}_world0_fine')
            d1=evaluate(1,probe,True,f'probe{index}_world1_fine')
        distance=float(np.linalg.norm(d1-d0)/.01)
        tasks.append(dict(probe=index,center=probe[0],radius=probe[1],contrast=probe[2],distance=distance))
        arrays.append((d0,d1));write(out/'ALL_TASKS.json',tasks)
        np.savez_compressed(out/f'probe{index}_PAIR.npz',data0=d0,data1=d1)
    numerical=all(max(r['world0_grid_distance'],r['world1_grid_distance'])<=.1 for r in rows)
    dcal=rows[0]['fine_distance'];dj=rows[1]['fine_distance'];dt=max(r['distance'] for r in tasks)
    gate=dict(identity='ORACLE_CONSTRUCTED_PHYSICAL_WITNESS',numerical_status='PASS_EMPIRICAL_GRID_CHECK' if numerical else 'NUMERICALLY_UNQUALIFIED',
        d_cal=dcal,d_joint_shifted=dj if kind=='joint' else None,d_T_sup_24=dt,
        task_argmax=max(tasks,key=lambda r:r['distance'])['probe'],
        numerical_pass=numerical,cal_physics_threshold=bool(dcal<=1 and dt>=3),
        joint_physics_threshold=bool(kind=='joint' and dj<=1),target_shift=.1 if kind=='joint' else 0,
        mass_relative_change=0 if kind=='joint' else None,
        original_minimum_feature_status='UNDEFINED_ORIGINAL_SCALAR; see independent geometry audit',
        strict_original_prior_GO=False,conditional_fixed_support_geometry=True,
        likelihood_not_posterior_mass=True,posterior_probability=None,totals=book.totals())
    write(out/'RESULTS.json',gate);print(json.dumps(gate),flush=True)

if __name__=='__main__':main(sys.argv[1])
