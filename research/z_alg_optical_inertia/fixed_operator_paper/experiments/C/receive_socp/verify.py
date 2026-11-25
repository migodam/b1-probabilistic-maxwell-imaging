"""Independently replay old conditional spectra and one-column supports at 192 bits."""
from collections import Counter
import csv
import sys
from time import perf_counter, process_time

from flint import arb, ctx
import numpy as np
from threadpoolctl import threadpool_limits

from solve import HERE, PAPER, SOURCE, bundle
from model import OLD, ball, hash_file, high, low, minor, norm, spectrum, star, write_json, read_json

sys.path.insert(0,str(PAPER/'protocol'))
from api import clique_cover_bound


def support(options,direction,idx,sigma):
    """Direct scalar square-root supports; no optimizer or T-factor is consulted."""
    w=ball(direction)
    size=norm(w)
    projected=[]
    radii=[]
    receipt=[]
    for e in options:
        center=minor(e.center,range(w.nrows()),[idx])
        beta=e.gram[idx,idx].real.nonnegative_part()
        r2=(star(w)*e.receive.solve(w))[0,0].real.nonnegative_part()
        spread=(beta*r2).nonnegative_part().sqrt()
        coord=(star(w)*center)[0,0].real
        pad=spread+sigma*size
        projected.append(coord)
        radii.append(spread)
        receipt.append(dict(clean_center_interval=[low(coord),high(coord)],
            clean_halfwidth_interval=[max(0.,low(spread)),high(spread)],
            noisy_support_interval=[low(coord-pad),high(coord+pad)]))
    margin=projected[1]-projected[0]-sum(radii)-2*sigma*size
    gap=low(arb(receipt[1]['noisy_support_interval'][0])-arb(receipt[0]['noisy_support_interval'][1]))
    return dict(checked=True,strict_edge=bool(margin>0 and gap>0),margin_interval=[low(margin),high(margin)],
                normalized_margin_interval=[low(margin/size),high(margin/size)] if size>arb('1e-8') else None,
                noisy_disjoint_interval_gap_lower=gap,per_class_support=receipt,
                direction_norm_interval=[max(0.,low(size)),high(size)],strict_unit_ball=bool(size<=arb(1)))


def labelled_control():
    labels=('A','B','C','D','E')
    edges={pair:dict(strict_edge=True,noisy_separation_gap_lower=.125) for pair in (
        ('A','B'),('A','C'),('B','C'),('D','E'))}
    edges['C','D']=dict(strict_edge=False,noisy_separation_gap_lower=.5)
    edges['A','E']=dict(strict_edge=True,noisy_separation_gap_lower=-.5)
    control=clique_cover_bound(labels,edges,dict(zip(labels,(1.,2.,3.,4.,5.),strict=True)))
    assert control['cliques']==[['A','B','C'],['D','E']]
    assert control['verified_edges']==4 and control['worst_survivors_upper']==2
    assert 8.<=control['worst_disjoint_area_upper']<8.00001
    return dict(labelled_only=True,invalid_edges_ignored=2,**control)


def main():
    wall,cpu=perf_counter(),process_time()
    pre=read_json(HERE/'proposals.json')
    for path,identity in pre['source_input_identities'].items():
        if hash_file(SOURCE/path)!=identity:
            raise ArithmeticError('Frozen old-only source input changed after proposal.')
    lookup={(x['scene'],x['contract'],x['source_index']):x for x in pre['rows']}
    rows=[]
    old_proofs=[]
    with np.load(HERE/'directions.npz',allow_pickle=False) as archive:
        for plan,g,nr,options,proofs in bundle():
            scene,contract=plan['scene'],plan['energy_contract']
            old_proofs.append(dict(scene=scene,contract=contract,verified_envelopes=proofs))
            for section in plan['batches']['1']['all_action_menu']:
                idx=section[0]
                row=lookup[scene,contract,idx]
                V=minor(g.ex_b,range(g.g.shape[0]),(*OLD,idx))
                source_floor=low(spectrum(star(V)*V)[0][0])
                if not source_floor>0.:
                    raise ArithmeticError('Old and new source directions are not strictly independent.')
                variants={}
                for kind in ('center','socp'):
                    key=row['key']+'__'+kind
                    variants[kind]=support(options,archive[key],idx,arb(nr)) if key in archive else dict(checked=False,strict_edge=False)
                c,o=variants['center'],variants['socp']
                if o['checked'] and not o['strict_unit_ball']:
                    raise ArithmeticError('Proposed SOCP direction is outside its verified norm constraint.')
                cm,sm=c['normalized_margin_interval'],o.get('normalized_margin_interval')
                difference=low(arb(sm[0])-arb(cm[1])) if cm and sm else None
                rows.append(dict(**row,source_gram_floor=source_floor,checks=variants,
                                 normalized_margin_improvement_lower=difference))
    def counts(subset):
        return dict(menu_actions=len(subset),center_checked=sum(x['checks']['center']['checked'] for x in subset),
            socp_checked=sum(x['checks']['socp']['checked'] for x in subset),
            center_strict=sum(x['checks']['center']['strict_edge'] for x in subset),
            socp_strict=sum(x['checks']['socp']['strict_edge'] for x in subset),
            new_strict=sum(x['checks']['socp']['strict_edge'] and not x['checks']['center']['strict_edge'] for x in subset),
            lost_strict=sum(not x['checks']['socp']['strict_edge'] and x['checks']['center']['strict_edge'] for x in subset),
            normalized_margin_strictly_improved_on_shared_edges=sum(
                x['checks']['center']['strict_edge'] and x['checks']['socp']['strict_edge'] and
                (x['normalized_margin_improvement_lower'] or 0.)>0. for x in subset))
    overall=counts(rows)
    assert overall['menu_actions']==120
    result=dict(status='COMPLETE_POST_STUDY_OLD_ONLY_COUNTERFACTUAL',source_round='diagnostic',
                independent_physical_objects=6,paired_contracts=2,angles_per_contract_condition=10,
                source_exclusions_from_existing_study_unchanged=True,future_data_read=False,verification_bits=ctx.prec,
                overall=overall,by_contract={s:counts([x for x in rows if x['contract']==s]) for s in ('added','fixed_total')},
                original_frozen_center_listed_statuses=dict(Counter(x['original_archive_status'] for x in rows)),
                solver_statuses=dict(Counter(x['proposal']['status'] for x in rows)),
                labelled_multiclass_control=labelled_control(),condition_prerequisites=old_proofs,rows=rows,
                proposal_wall_seconds=pre['wall_seconds'],proposal_cpu_seconds=pre['cpu_seconds'],
                verification_wall_seconds=perf_counter()-wall,verification_cpu_seconds=process_time()-cpu)
    write_json(HERE/'results.json',result)
    with (HERE/'results.csv').open('w',newline='') as f:
        out=csv.writer(f)
        out.writerow(('scene','contract','source_angle','old_archive_status','center_checked','center_strict','center_margin_low',
                      'socp_checked','socp_strict','socp_margin_low','normalized_margin_improvement_lower'))
        for x in rows:
            c,o=[x['checks'][k] for k in ('center','socp')]
            out.writerow((x['scene'],x['contract'],x['source_angle'],x['original_archive_status'],c['checked'],c['strict_edge'],
                          c['margin_interval'][0],o['checked'],o['strict_edge'],o.get('margin_interval',[None])[0],
                          x['normalized_margin_improvement_lower']))
    print(result['overall'])
    print(result['by_contract'])


if __name__=='__main__':
    ctx.prec=192
    ctx.threads=1
    with threadpool_limits(1):
        main()
