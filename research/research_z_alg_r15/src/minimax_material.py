"""Two-point minimax material error from the physical MIMO complex mean difference."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

OUT=Path(__file__).resolve().parents[1]/'results'


def conditional_material_risk(gap, distance, pe):
    # Along the real whitened mean direction, the binary populations are N(±b/2,1).
    # Under the + hypothesis, log posterior odds for + are b^2/2 + b Z.
    b=np.sqrt(2.)*gap
    def loss(z):
        logit=b*(z+b/2)
        return np.exp(-.5*z*z)/np.sqrt(2*np.pi)*expit(logit)*expit(-logit)
    value,_=quad(loss,-45.,45.,epsabs=1.e-28,epsrel=1.e-11,points=[-b/2,0.],limit=130)
    return dict(optimal_binary_error=pe,lower_expected_material_rmse=distance*pe,
                lower_root_expected_material_mse=distance*np.sqrt(value))


def main():
    physical=json.loads((OUT/'witness_transfer.json').read_text())['cases']
    maps=json.loads((OUT/'material_weak_modes.json').read_text())['cases']
    result=[]
    for p in physical:
        m=next(x for x in maps if x['hard_design']==p['hard_design'] and x['easy_design']==p['easy_design'])
        dif=np.array(m['plus'])-np.array(m['minus'])
        distances={'full_576':float(np.linalg.norm(dif)/24.),'patch_88':float(np.linalg.norm(dif)/np.sqrt(len(dif)))}
        rows={}
        for side in ('hard','easy'):
            v=p['grid']['72'][side];rows[p[side+'_design']]={area:conditional_material_risk(v['gap_over_sigma'],d,v['equal_prior_error']) for area,d in distances.items()}
        result.append(dict(material_pair_hard_for=p['hard_design'],material_rms_distances=distances,designs=rows))
    (OUT/'minimax_material.json').write_text(json.dumps(dict(cases=result,assumptions='Finite 72 VIE means and known gain, independent circular Gaussian receiver noise, fixed acquisition, and any estimator on a parameter set containing these endpoints. The numerical means have no outward enclosure.'),indent=2)+'\n')


if __name__=='__main__':
    main()
