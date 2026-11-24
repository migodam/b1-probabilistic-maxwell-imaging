"""Offline truth-in-box and no-false-prune ledger, after the separate online seal."""

import argparse
import json
import sys

import numpy as np
from threadpoolctl import threadpool_limits

from branches import as_branch,central_safety_control,expand
from holdout_public import independent_screen
from native_box import HERE,NativeBox,clock,elapsed

sys.path.insert(0,str(HERE.parent))
from holdout.data import APERTURES,DEFAULT_ROOT,get_native_case,get_native_model,planned_cases,write_json


def oracle_support(truth):
    occupied=np.any((truth.reshape(96,96)>0).reshape(8,12,8,12),axis=(1,3))
    return as_branch('offline_true_tile_cover',np.where(occupied,max(4.,float(truth.max())),0.))


def row(case,name,branch,material,reference_feasible,result):
    caps=expand(branch,96)
    truth_in_box=bool(np.all((material>=0)&(material<=caps)))
    false=bool(reference_feasible and truth_in_box and result['decision']=='PRUNE')
    return dict(case=case,arm=name,truth_in_box=truth_in_box,native_truth_ball_verified=reference_feasible,
                counted_true_branch=truth_in_box and reference_feasible,false_prune=false,
                truth_peak=float(material.max()),tile_caps=branch.caps,result=result)


def tally(rows):
    aggregate={}
    for item in rows:
        key='/'.join((item['arm'],item['result']['branch']))
        values=aggregate.setdefault(key,dict(rows=0,true_covered=0,false_pruned=0,actually_pruned=0,
                                             materials=set(),covered_materials=set(),false_materials=set()))
        values['rows']+=1
        values['true_covered']+=item['counted_true_branch']
        values['false_pruned']+=item['false_prune']
        values['actually_pruned']+=item['result']['decision']=='PRUNE'
        identity=item['case'].split('__',1)[0]
        values['materials'].add(identity)
        if item['counted_true_branch']:
            values['covered_materials'].add(identity)
        if item['false_prune']:
            values['false_materials'].add(identity)
    for values in aggregate.values():
        for key in ('materials','covered_materials','false_materials'):
            values['independent_'+key]=len(values.pop(key))
    return aggregate


def main(args):
    if not args.online_seal.is_file():
        raise FileNotFoundError('Require the online cohort seal before scoring protected references.')
    controls=json.loads(args.broad.read_text())['broad']['cases']
    score=json.loads(args.truth_audit.read_text())
    truths={r['observation_id']:bool(r['true_material_in_padded_ball'])
            for r in score['observations']}
    cases=list(planned_cases())
    if set(truths)!={r['case'] for r in cases}:
        raise ValueError('Require the full paired native truth-ball audit before safety scoring.')
    if any(not (args.predictions/(entry['case']+'__pruned.json')).is_file() for entry in cases):
        raise ValueError('Expected a frozen pruned receipt for every observation.')
    start=clock()
    m_cache={}
    rows=[]
    for entry in cases:
        name,oid=entry['case'],entry['object_id']
        if oid not in m_cache:
            with np.load(args.root/'offline/objects'/f'{oid}.npz',allow_pickle=False) as labels:
                m_cache[oid]=labels['m_native'].reshape(-1)
        ref=truths.get(name,False)
        m=m_cache[oid]
        online=json.loads((args.predictions/(name+'__pruned.json')).read_text())
        for data in online['certificates']:
            if 'tile_caps' not in data:
                continue
            b=as_branch(data['branch'],data['tile_caps'])
            rows.append(row(name,'online_primary',b,m,ref,data))
        rows.append(row(name,'predeclared_central',central_safety_control(),m,ref,controls[name]))
    summary=lambda complete:dict(online_seal=str(args.online_seal),truth_audit=str(args.truth_audit),
                                  oracle_checks_complete=complete,summary=tally(rows),rows=rows,**elapsed(start))
    write_json(args.output,summary(False))
    print('actual online and predeclared central coverage',tally(rows),flush=True)
    if args.coverage_only:
        return
    g={}
    for label in APERTURES:
        p=get_native_model(label,args.root)
        g[label]=NativeBox(p['spectrum'],p['S'],p['E'])
    cache={}
    for entry in cases:
        name,oid=entry['case'],entry['object_id']
        m=m_cache[oid]
        ref=truths.get(name,False)
        branch=oracle_support(m)
        case,p=get_native_case(name,args.root)
        geom=g[case['array']]
        radius=float(p['delta'])
        key=case['array'],branch.caps
        if key in cache:
            evaluation=dict(branch=branch.name,decision='RETAIN',residual_lower=0.,radius=radius,
                            reused_rank_exhaustion_from=cache[key])
        else:
            evaluation=geom.score(p['Y'],radius,branch,max_proofs=1,max_rank=850)
            if independent_screen(evaluation,geom.ns):
                cache[key]=name
        rows.append(row(name,'offline_oracle',branch,m,ref,evaluation))
        write_json(args.output,summary(False))
    write_json(args.output,summary(True))
    print('complete safety',tally(rows),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=type(HERE),default=DEFAULT_ROOT)
    parser.add_argument('--predictions',type=type(HERE),default=HERE.parent/'benchmark/predictions')
    parser.add_argument('--online-seal',type=type(HERE),required=True)
    parser.add_argument('--truth-audit',type=type(HERE),required=True)
    parser.add_argument('--broad',type=type(HERE),default=HERE/'holdout_public_broad.json')
    parser.add_argument('--output',type=type(HERE),default=HERE/'holdout_offline_safety.json')
    parser.add_argument('--coverage-only',action='store_true')
    with threadpool_limits(1):
        main(parser.parse_args())
