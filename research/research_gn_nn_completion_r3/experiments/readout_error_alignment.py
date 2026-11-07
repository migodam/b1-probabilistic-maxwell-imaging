"""Post-seal diagnostic of the conditional weighted readout normal."""
from pathlib import Path
import argparse,json
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    scores=json.loads(a.scores.read_text())
    seal=json.loads(Path(scores['seal_path']).read_text())
    if seal['status']!='ALL_PLANNED_PREDICTIONS_SEALED':
        raise ValueError('Only frozen offline diagnostics')
    study=Path(seal['study']);data=Path(seal['data_root'])
    first={}
    for r in scores['rows']:first.setdefault(r['obs_id'],r)
    rows=[]
    for obs,r in first.items():
        with np.load(study/'anchors'/f'{obs}.npz')as z:
            m0=z['m0'];q=z['Q0']
        with np.load(data/r['truth_path'])as z:err=z['m']-m0
        d=np.sum(abs(q)**2,axis=1)
        fraction=float((err@d)**2/max((err@err)*(d@d),1e-300))
        rows.append(dict(obs_id=obs,parent_id=r['parent_id'],level=r.get('level','low'),
                         noise=r['noise'],distribution=r.get('distribution','ID'),
                         GN_GT_material_error_energy_fraction_along_weighted_readout_normal=fraction))
    groups=[]
    for level in ('low','stress'):
        for noise in (0.,.001):
            v=[r['GN_GT_material_error_energy_fraction_along_weighted_readout_normal']
               for r in rows if r['level']==level and r['noise']==noise]
            if v:groups.append(dict(level=level,noise=noise,n=len(v),mean_fraction=float(np.mean(v)),
                                    median_fraction=float(np.median(v)),max_fraction=float(max(v))))
    a.output.write_text(json.dumps(dict(scope='OFFLINE after TEST predictions sealed; no fit, basis selection, new LU or RHS',
        interpretation='d_i=||Q0_i||² is a conditional readout-annihilator direction only under the stated fixed-operator inclusion conditions; this is not the physical material Jacobian nullspace, nor proof of global range.',groups=groups,rows=rows),indent=2))
    print(json.dumps(groups))

if __name__=='__main__':main()
