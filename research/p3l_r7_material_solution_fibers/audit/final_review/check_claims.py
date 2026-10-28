"""Final claim readback: sealed public predictions/receipts only, zero GT/physics."""
import json
from pathlib import Path
import numpy as np
from scipy.special import ndtr
from scipy.stats import chi2

P=Path(__file__).resolve().parents[2]
sigma=.00119211798514519
def white(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()*np.sqrt(2)/sigma
def load(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def energy(z,y):return float(np.linalg.norm(white(z-y))**2/2)
summary={a:dict(executed=0,competitive=0,competitive_pair_le1=0) for a in ('A_FROZEN','B_FULLJ','C_JETS','D_RANDOM')}
pairs=[];first=[];paths=[];cap=[]
for i in range(8):
    folder=P/'E2/parents'/f'R7_LOCKED_{i:03d}'
    base=load(folder/'ANCHOR.npz');y=load(folder/'OBSERVATION.npz')['observation']
    geometry=load(folder/'BASE_GEOMETRY.npz');j=load(folder/'COMPENSATED_JETS.npz')
    c1=j['c1'];A=j['A']
    first.append(dict(parent=folder.name,norm=float(np.linalg.norm(c1)),
        strong_projectable_fraction=float(np.linalg.norm(A@np.linalg.lstsq(A,c1,rcond=1e-10)[0])/np.linalg.norm(c1)),
        chain_rule_relative=float(np.linalg.norm(c1-base['J']@j['v'])/np.linalg.norm(c1)),
        material_speed=float(np.linalg.norm(np.r_[j['m1'].real,j['m1'].imag])/24)))
    status=json.loads((folder/'ANCHOR_STATUS.json').read_text())
    cap.append(dict(parent=folder.name,fresh=status['fresh_state_attempts'],optimizer=status['optimizer']))
    assert status['fresh_state_attempts']<=64
    for arm,counts in summary.items():
        path=folder/arm/'RMS_0.0200/ENDPOINT.npz'
        if not path.exists():continue
        p=load(path);counts['executed']+=1
        material=p['material'] if 'material' in p else p['chi']
        base_material=base['material'] if 'material' in base else base['chi']
        dm=np.r_[(material-base_material).real,(material-base_material).imag]/24
        unit=geometry['random_material'] if arm=='D_RANDOM' else geometry['weak_material']
        phi=energy(p['prediction'],y)
        d=float(np.linalg.norm(white(p['prediction']-base['prediction'])))
        legal=np.isfinite(p['q']).all() and (p['q'][64:]>=0).all() and (p['q'][64:]<=1).all()
        good=(abs(np.linalg.norm(dm)-.02)<=.0004 and unit@dm>=.01 and legal
              and phi<=chi2.ppf(.95,512)/2 and phi-energy(base['prediction'],y)<=2)
        counts['competitive']+=int(good)
        counts['competitive_pair_le1']+=int(good and d<=1)
        if arm=='B_FULLJ':pairs.append(dict(parent=folder.name,pair_L2=d,competitive=bool(good)))
    for name in ('B_FULLJ_CHORD','C_RAW_CUBIC'):
        file=folder/'paths'/name/'PATH.json'
        if not file.exists():continue
        doc=json.loads(file.read_text());nodes=doc['nodes']
        row=dict(parent=folder.name,path=name,nodes=len(nodes),qualified=doc['sampled_path_qualified'])
        assert len(nodes)==5
        if name=='C_RAW_CUBIC':
            t1=nodes[-1];node=load(file.parent/'NODE_1.00.npz');h=t1['h']
            actual=white(node['prediction']-base['prediction'])
            predicted=h*j['c1']+h*h*j['c2']+h**3*j['c3']
            row.update(endpoint_prediction_relative=float(np.linalg.norm(actual-predicted)/np.linalg.norm(actual)),
                       endpoint_seam_crossings=t1['seam_crossings'])
        paths.append(row)

resource=json.loads((P/'RESOURCE_SUMMARY.json').read_text())
ledger_rows=[json.loads(line) for line in (P/'RESOURCE_LEDGER.jsonl').read_text().splitlines()]
charges={k:sum(row[k] for row in ledger_rows) for k in ('RHS','states','LU')}
assert charges==dict(RHS=37456,states=1325,LU=1325)
assert len([r for r in ledger_rows if r['event']=='JOB_BEGIN'])==19
assert sum(resource['categories_RHS'].values())==37456
assert resource['forward_RHS']+resource['derivative_RHS']==37456
assert charges['RHS']==resource['new_complex_RHS']
assert json.loads((P/'CLOSURE.json').read_text())['total_RHS']==37456

analysis=json.loads((P/'analysis/MATERIAL_PHYSICS.json').read_text())
risk_errors=[]
for r in analysis['two_world_tasks']:
    d=r['whitened_pair_L2'];tv=2*ndtr(d/2)-1
    lower=np.asarray(r['tasks_difference'])**2*(1-tv)/8
    risk_errors.append(max(abs(d*d/2-r['two_fixed_Gaussian_worlds_KL']),
        abs(tv-r['two_fixed_Gaussian_worlds_TV']),
        float(np.max(abs(lower-r['scalar_two_point_minimax_lowerbounds'])))))
assert max(risk_errors)<1e-12
report=dict(status='CLAIM_NUMBERS_PASS',scope='sealed endpoint/current-jet arrays and receipts only; no GT, no physical imports',
    summary_002=summary,B_pairs=pairs,first_order=first,paths=paths,anchors=cap,
    journal_charges=charges,cost_categories=resource['categories_RHS'],
    scalar_fixed_pair_formula_max_error=max(risk_errors),new_Maxwell_RHS=0,GT_access=False)
(P/'audit/final_review/CLAIM_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],summary=summary,charges=charges,
    c1_norm_range=[min(x['norm'] for x in first),max(x['norm'] for x in first)],
    first_order_strong_projectable_fraction_range=[min(x['strong_projectable_fraction'] for x in first),max(x['strong_projectable_fraction'] for x in first)],
    scalar_fixed_pair_formula_max_error=max(risk_errors)),indent=2))
