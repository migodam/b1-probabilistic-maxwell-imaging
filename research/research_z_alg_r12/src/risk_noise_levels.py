"""Fresh noise draws test the physical bandwidth change when the noise level changes."""
import json
import numpy as np
from freq_fold import candidate_basis
from freq_model import FrequencyTM
from mirror_spectrum import Mirror,CENTER,V,STRONG,ROOT
from risk_design import study

if __name__=='__main__':
    old=json.loads((ROOT/'results/risk_design.json').read_text())
    rows=[]
    for level_index,level in enumerate([.15,.90]):
        for geom_index,r in enumerate(old):
            pair=np.asarray(r['pair']);a=r['amplitude']
            if r['label']=='mirror':
                model=Mirror([1,2]);basis=(CENTER,V,STRONG,a,True)
            else:
                model=FrequencyTM(r['geometry']);basis=candidate_basis(model,pair)
            rows.append(study(r['label'],model,basis,pair,24,308737+200*geom_index+5000*level_index,level))
            (ROOT/'results/risk_noise_levels.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
