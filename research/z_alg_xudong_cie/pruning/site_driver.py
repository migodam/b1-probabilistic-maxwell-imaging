"""Blindly evaluate all 16 exact-air tile hypotheses on compact native site controls."""

import argparse
from datetime import datetime,timezone
import json

import numpy as np
from threadpoolctl import threadpool_limits

from native_box import HERE,clock,elapsed
from screen import Checks
from site_menu import SiteMenu

DEFAULT_PUBLIC=HERE.parents[2]/'.artifacts/z_alg_xudong_cie/site_control/public'


def write(path,data):
    tmp=path.with_suffix('.pending')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)


def run(args):
    start=clock()
    manifest_path=args.public/'manifest.json'
    manifest=json.loads(manifest_path.read_text())
    requested=manifest['cases']
    if len(requested)!=16 or len({entry['object_id'] for entry in requested})!=4:
        raise ValueError('Require four independent objects with exactly sixteen paired views.')
    geometry={}
    prep={}
    with np.load(args.public/'geometry96.npz',allow_pickle=False) as archive:
        for ap in ('PW180_10x20','PW360_10x20'):
            t=clock()
            checks=Checks(archive['spectrum'],archive[ap+'__S'],archive[ap+'__B'])
            if checks.aperture!=ap:
                raise ValueError('Site geometry differs from the frozen native source/receiver arrays.')
            geometry[ap]=SiteMenu(checks)
            prep[ap]=dict(**elapsed(t),ordinary_born_setup=geometry[ap].born_setup,
                          inherited_geometry_setup=checks.setup_cost,
                          translated_tile_gap_proof='charged on the first actual tile call')
    out=dict(scope='public compact-site control; prior to any generating-location read',
             manifest=str(manifest_path),predictions={},aperture_preparations=prep,
             model='All sixteen 4x4 tiles, inside contrast <= 4 and exactly zero exterior, with the frozen native Checks API.',
             classifier='No forced winner. Every tile with no strict exclusion is retained.',
             born_metric='Data residual against first-Born unit tile, with one optimum complex scale per tile.',
             state='IN_PROGRESS')
    for entry in requested:
        t=clock()
        source=(args.public/entry['measurement_file']).resolve()
        if not source.is_relative_to(args.public.resolve()):
            raise ValueError('Measurement path must stay inside the public root.')
        with np.load(source,allow_pickle=False) as observed:
            if set(observed.files)!={'Y','delta','raw_radius'}:
                raise ValueError('The online site model expects only measurement and two noise radii.')
            y,radius=observed['Y'],float(observed['delta'])
        result=geometry[entry['array']].score(y,radius)
        out['predictions'][entry['case']]=dict(object_id=entry['object_id'],aperture=entry['array'],
                                             noise_relative=entry['noise_relative'],
                                             source=str(source.relative_to(args.public.resolve())),
                                             **result,full_case_time=elapsed(t))
        out['online_cost_to_date']=elapsed(start)
        write(args.output,out)
        print(entry['case'],'survivors',result['survivors'],'one-sided',len(result['receive_only_survivors']),
              len(result['illumination_only_survivors']),'Born best',result['born_ranking'][0]['tile'],flush=True)
    out['state']='SEALED_PUBLIC_SITE_OUTPUTS'
    out['sealed_at']=datetime.now(timezone.utc).isoformat()
    out['online_cost_to_date']=elapsed(start)
    write(args.output,out)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public',type=type(HERE),default=DEFAULT_PUBLIC)
    parser.add_argument('--output',type=type(HERE),default=HERE/'site_online.json')
    with threadpool_limits(1):
        run(parser.parse_args())
