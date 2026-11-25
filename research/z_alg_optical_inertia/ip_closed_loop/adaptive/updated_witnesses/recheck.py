"""Replay public old and fresh noise balls at a location's physical witness."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

OUT = Path(__file__).resolve().parent
ADAPTIVE = OUT.parent
MAPS = ADAPTIVE.parent/'maps'
sys.path.insert(0,str(MAPS))
from common import FREEZE, Geometry, clock, elapsed, geom_path, pattern, place, save_json
from sector_map import SectorGeometry, matrix_epsilon
from certified_map import fnorm
from witness import certify_witness, forecast, norm_sensor

BACKGROUND = FREEZE['model']['background_upper']
CORE_LO,CORE_HI = FREEZE['material_classes']['primary_core_range']


def read(path):
    return json.loads(path.read_text())


def verify_material(g,x,cells,profile):
    upper=np.full(g.n,BACKGROUND)
    lower=np.zeros(g.n)
    lower[cells]=CORE_LO*profile
    upper[cells]=CORE_HI*profile
    if not np.all(np.isfinite(x)) or not np.all((x>=lower)&(x<=upper)):
        raise ValueError('Witness exceeds the registered material bounds.')


def recheck_using_inverse(sector,x,incident,current,data,radius,inverse):
    """Reuse the same-material source-independent I-GX inverse from certify_witness."""
    t=clock()
    g=sector.g
    wave=sector.apply(current,g.symbol)
    defect=current-x[:,None]*(incident+wave)
    xmax=float(np.max(x))
    error_gj=sector.rel_fft_bound*sector.full_g*fnorm(current)
    arithmetic=matrix_epsilon(g.n)*(1+fnorm(current)+xmax*(fnorm(incident)+fnorm(wave)))
    defect_bound=fnorm(defect)+xmax*error_gj+arithmetic
    current_bound=(1+xmax*inverse*sector.full_g)*defect_bound
    observed=g.S@current
    roundoff=matrix_epsilon(g.n)*(1+fnorm(data)+fnorm(g.S)*fnorm(current))
    residual=fnorm(data-observed)+norm_sensor(sector)*current_bound+roundoff
    return dict(residual_upper=float(residual),noise_radius=float(radius),
                residual_ratio_upper=float(residual/radius),physical_state_defect_upper=float(defect_bound),
                verified=bool(residual<radius),**elapsed(t))


def make_report(side, do_fit):
    start=clock()
    g=Geometry(geom_path(side))
    sector=SectorGeometry(g,BACKGROUND)
    g_base=g.B.copy()
    stock=read(MAPS/'public'/'acquisition.json')['observations']
    receipts={}
    rows={}
    for path in (MAPS/'results').glob(f'witnesses_{side}_*.json'):
        for oid,v in read(path)['witnesses'].items():
            if '__weak__' in oid:
                receipts[oid]=v.get(str(BACKGROUND),{}).get('witness')
    for oid,info in stock.items():
        if info['side']!=side or info['scene']!='weak':
            continue
        folder=ADAPTIVE/oid
        plan=read(folder/'proposal.json')
        old=receipts[oid]
        if old is None:
            raise ValueError(f'Missing archived physical witness for {oid}.')
        row_tick=clock()
        old_data=np.load(MAPS/'public'/f'{oid}.npz',allow_pickle=False)
        Yold,r=old_data['Y'],float(old_data['r'])
        if float(old['noise_radius'])!=r:
            raise ArithmeticError('Archived witness radius does not match public old observation.')
        source=np.load(MAPS/'results'/old['file'],allow_pickle=False)
        x,cells,profile,current=[source[z] for z in ('material','active_cells','geometric_profile','current')]
        verify_material(g,x,cells,profile)
        coordinate=[int(c) for c in source['centre']]
        center=coordinate[1]*side+coordinate[0]
        if coordinate!=old['centre']:
            raise ArithmeticError('Stored witness coordinate differs from archived record.')
        expected_cells,expected_profile=place(side,coordinate,*pattern(side,info['shape']))
        if not np.array_equal(cells,expected_cells) or not np.array_equal(profile,expected_profile):
            raise ArithmeticError('Stored witness footprint differs from the registered location class.')
        old_unchanged=recheck_using_inverse(sector,x,g_base,current,Yold,r,old['full_field_inverse_upper'])
        if not old_unchanged['verified']:
            raise ArithmeticError('Original archived old material did not replay.')
        row=dict(centre=coordinate,archived_witness=old['file'],stopped=plan['stopped'],
                 old_original=old_unchanged,arms={})
        if plan['stopped']:
            if center not in plan['survivors']:
                raise ArithmeticError('Stopped centre did not retain its original witness.')
            row.update(**elapsed(row_tick))
            rows[oid]=row
            continue
        proposed=np.load(folder/'proposal.npz',allow_pickle=False)
        fresh=np.load(folder/'future.npz',allow_pickle=False)['Y']
        update=np.load(folder/'update.npz',allow_pickle=False)
        survivors=update['ids']
        margins=update['margins']
        if not np.isclose(plan['old_radius'],r,rtol=0,atol=0):
            raise ArithmeticError('Old radius changed between old and fresh stages.')
        new_radius=plan['new_radius']
        g.B=proposed['B']
        measurement,current,code=g.forward(x)
        if code:
            raise ArithmeticError(f'Forward replay failed for {oid}.')
        for j,name in enumerate(plan['policies']):
            tick=clock()
            Ynew=fresh[:,j,None]
            incident=proposed['B'][:,j,None]
            point=current[:,j,None]
            fresh_check=recheck_using_inverse(sector,x,incident,point,Ynew,new_radius,old['full_field_inverse_upper'])
            centre_was_retained=bool(center in survivors[~(margins[:,j]>0)])
            checks=dict(original_material_new=fresh_check,returned_centre_survives=centre_was_retained,
                        original_verified_both=fresh_check['verified'])
            if fresh_check['verified'] and not centre_was_retained:
                raise ArithmeticError('Updated map contradicts its original-operator physical witness.')
            if do_fit and not fresh_check['verified'] and centre_was_retained:
                checks['fitted']=fit_once(g,sector,oid,name,x,cells,profile,coordinate,g_base,incident,Yold,r,Ynew,new_radius)
            checks.update(**elapsed(tick))
            row['arms'][name]=checks
        row.update(**elapsed(row_tick))
        rows[oid]=row
        print(oid.split('__')[1],[(a,round(v['original_material_new']['residual_ratio_upper'],4),
                                    v.get('fitted',{}).get('joint_verified'))
                                   for a,v in row['arms'].items()],flush=True)
    g.B=g_base
    summary=dict(side=side,readings='public old and acquired new complex data; no generator arrays read',
                 unchanged_material_inverse='archived original-operator bound for I-GX, same model and material; source-independent',
                 rows=rows,**elapsed(start))
    save_json(OUT/f"{'audit' if do_fit else 'unchanged'}_{side}.json",summary)


def fit_once(g,sector,oid,name,x,cells,profile,coordinate,Bold,Bnew,Yold,r,Ynew,new_radius):
    t=clock()
    result=dict(joint_verified=False)
    original_B=g.B
    ref=float(np.median(np.delete(x,cells)))
    old_count=Yold.shape[1]
    try:
        for preference in (1.,1.5,2.):
            ratio=r*preference/new_radius
            g.B=np.column_stack((Bold,ratio*Bnew))
            joined=np.column_stack((Yold,ratio*Ynew))
            tentative=forecast(g,x,cells,profile,joined,r,BACKGROUND)
            if tentative is None:
                continue
            candidate,current,_=tentative
            verify_material(g,candidate,cells,profile)
            old_curr=current[:,:old_count]
            new_curr=current[:,old_count:]/ratio
            num_old=fnorm(Yold-g.S@old_curr)/r
            num_new=fnorm(Ynew-g.S@new_curr)/new_radius
            result.setdefault('trials',[]).append(dict(relative_new_weight=preference,old_numeric=float(num_old),new_numeric=float(num_new)))
            if max(num_old,num_new)>=1:
                continue
            g.B=Bold
            audit=certify_witness(sector,candidate,cells,old_curr,Yold,r,ref)
            if audit is None:
                continue
            fresh=recheck_using_inverse(sector,candidate,Bnew,new_curr,Ynew,new_radius,audit['full_field_inverse_upper'])
            if not fresh['verified']:
                continue
            folder=OUT/'materials'
            folder.mkdir(exist_ok=True)
            target=folder/f'{oid}__{name}.npz'
            np.savez_compressed(target,material=candidate,centre=coordinate,active_cells=cells,
                                old_current=old_curr,new_current=new_curr)
            result.update(joint_verified=True,old=audit,new=fresh,material=target.relative_to(OUT).as_posix())
            break
    finally:
        g.B=original_B
    return dict(result,**elapsed(t))


def summarize():
    batch=[read(OUT/f'audit_{side}.json') for side in (32,96)]
    acquired=0
    evaluated=0
    stopped=0
    unchanged=0
    refitted=0
    by_policy={}
    max_old=0.
    max_new=0.
    for result in batch:
        for row in result['rows'].values():
            if row['stopped']:
                stopped+=1
                continue
            acquired+=1
            for policy,check in row['arms'].items():
                evaluated+=1
                numbers=by_policy.setdefault(policy,dict(unchanged_verified=0,additional_verified_after_fit=0,
                                                          verified=0,retained_centres=0))
                numbers['retained_centres']+=int(check['returned_centre_survives'])
                fitted=check.get('fitted',{})
                reused=check['original_verified_both']
                recovered=fitted.get('joint_verified',False)
                unchanged+=int(reused)
                refitted+=int(recovered)
                numbers['unchanged_verified']+=int(reused)
                numbers['additional_verified_after_fit']+=int(recovered)
                numbers['verified']+=int(reused or recovered)
                if reused:
                    old=row['old_original']
                    new=check['original_material_new']
                elif recovered:
                    old=fitted['old']
                    new=fitted['new']
                else:
                    continue
                max_old=max(max_old,old['residual_ratio_upper'])
                max_new=max(max_new,new['residual_ratio_upper'])
    receipt=dict(total_maps=acquired+stopped,already_stopped_old_witnesses=stopped,
                 new_acquired_maps=acquired,paired_policy_observations=evaluated,
                 original_verified_joint=unchanged,updated_material_verified_joint=refitted,
                 joint_total=unchanged+refitted,per_policy=by_policy,
                 maximum_verified_old_radius_ratio=max_old,maximum_verified_new_radius_ratio=max_new,
                 final_audit_wall_seconds=sum(r['wall_s'] for r in batch),
                 final_audit_cpu_seconds=sum(r['cpu_s'] for r in batch),
                 scope='One material witness at the retained centre for each actual one-shot policy, verified separately in its old and new noise balls. Separate policies can have different refitted materials; no GT access, no certification of every remaining centre.')
    save_json(OUT/'summary.json',receipt)
    print(receipt,flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('side',type=int,choices=(32,96))
    parser.add_argument('--fit',action='store_true')
    args=parser.parse_args()
    make_report(args.side,args.fit)
    if args.fit and all((OUT/f'audit_{side}.json').exists() for side in (32,96)):
        summarize()


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
