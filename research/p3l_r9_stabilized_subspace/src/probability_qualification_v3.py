"""Zero-physics postprocessing producing canonical numeric qualification v3.

Preserves previous qualification files, fixes source-refinement inclusion and
labels finite-refinement denominator/omission as proxies rather than proof.
"""
import json
from pathlib import Path
P=Path(__file__).resolve().parents[1];p=P/'results/probability'
v=json.loads((p/'final_numerical_comparison.json').read_text());old=json.loads((p/'integration_results.json').read_text())
for name,r in v.items():
    src=old[name]['source_coordinate_integrals'];gh21=json.loads((p/f'{name}_reduced_GH21.json').read_text());gh31=r['reference'];zz=[x['scaled_evidence'] for x in [*src,gh21,gh31]]
    denominator=min(zz)-2*(max(zz)-min(zz));assert denominator>0
    q=r['reference_qualification'];q['scaled_evidence_numerical_denominator_proxy']=denominator;q.pop('scaled_evidence_numerical_lower_estimate',None)
    q['numerical_omitted_posterior_upper_proxy']=r['whole_support_audit']['numerical_scaled_omitted_evidence_upper']/denominator;q.pop('numerical_omitted_posterior_upper',None)
    q['numerical_denominator_formula']='min(all source/GH refinements) - 2*(max(all)-min(all)); positive conservative finite-refinement proxy, NOT a guaranteed lower bound'
    q['qualification']='NUMERICAL_REFERENCE: finite cross-coordinate/task refinement plus entire-support float majorant audit. Neither the denominator proxy, omission proxy nor reported cross-refinement uncertainty is a rigorous certificate.'
    q['all_relative_scaled_evidence_spread']=(max(zz)-min(zz))/gh31['scaled_evidence']
    for order,rr in r['two_factor_physical_reuse'].items():
        rr['shared_setup_plus_integral_wall']=rr.pop('shared_cold_numerical_path_wall')
        rr['interpretation']='Single-run same quadrature physical reuse components only. Setup plus integral does not include common fresh profile, qualification reference, all validation, file reads or process imports; it is NOT a full cold-start timing. Total campaign receipt includes qualification costs separately.'
    (p/f'{name}_reference_qualification_v3.json').write_text(json.dumps(q,indent=2)+'\n')
(p/'final_numerical_comparison_v3.json').write_text(json.dumps(v,indent=2)+'\n')
