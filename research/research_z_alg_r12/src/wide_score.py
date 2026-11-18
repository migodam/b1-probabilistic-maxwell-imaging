"""Material structure and cost of the saved, same-data wide-band baselines."""
from pathlib import Path
import json
import numpy as np
from high_score import ENDPOINTS, material_scores, truth_for

OUT=Path(__file__).resolve().parents[1]/'results/freq_high'

def main():
    table=[]
    for label in ('ring','shell','bridge'):
        truth,support,amplitude,_=truth_for(ENDPOINTS[label].case_id)
        for budget in (1440,7200):
            stem=OUT/f'{label}_wide_{budget}'
            records=json.loads(Path(str(stem)+'_solves.json').read_text())
            with np.load(str(stem)+'_fits.npz') as mats:
                for name,v in records.items():
                    elapsed=v.get('elapsed',sum(x['elapsed'] for x in v.get('stages',[])))
                    table.append(dict(case=label,method=name,budget=budget,seconds=elapsed,
                                      LU=v['counts']['LU'],RHS=v['counts']['actual_complex_RHS'],
                                      **material_scores(mats[name],truth,support,amplitude)))
    (OUT/'structural_scores.json').write_text(json.dumps(table,indent=2)+'\n')

if __name__=='__main__':
    main()
