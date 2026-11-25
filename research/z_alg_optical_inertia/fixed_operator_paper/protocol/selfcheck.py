"""Execute every protocol entry on a small control or a frozen experimental input."""
import argparse
import json
from pathlib import Path
import shutil
import tempfile
from time import perf_counter

from flint import ctx
import numpy as np
from threadpoolctl import threadpool_limits

import api


def execute():
    started=perf_counter()
    proof={}
    g=np.zeros((4,4))
    g[0,1]=g[1,0]=g[2,3]=g[3,2]=.125
    model=api.ModelContract(g,np.array([[1.,0,0,0]]),np.eye(4)[:,:2],(2,2),'synthetic check')
    box=api.MaterialNode('witness-box','control',np.zeros(4),np.full(4,.5))
    noise=api.NoiseContract(.0625,'known deterministic synthetic radius')
    old=api.Acquisition(model,[[1],[0]],[[2.]],noise,'unit receiver coordinates',(1.,),np.eye(2))
    other=api.Acquisition(model,[[0],[1]],[[0]],noise,'unit receiver coordinates',(1.,),np.eye(2))
    weights=(1.,0.)
    base=api.assemble_border(model,box,weights,old)
    assert base.K.nrows()==4 and base.L.ncols()==1 and base.C.nrows()==1
    proof['26']={'dimensions':[base.K.nrows(),base.C.nrows()],'hashes':base.hashes}

    rows=api.local_pair_rows(model,[0],old)
    T,E=rows
    assert api.trace.arr(T)[0,0]==1 and api.trace.arr(E)[1,0]==1
    proof['27']={'ordered_rows':T.nrows(),'source_affine_entry':1.}
    kernel=api.build_fixed_kernel(base,rows)
    assert kernel.reduced.eliminated==2
    proof['28']={'retained':kernel.reduced.dim,'eliminated_positive':kernel.reduced.eliminated,
                  'proof':kernel.reduced.gram_rank_scope}

    delta=np.diag([0.,.25])
    reduced=api.update_box(kernel,delta,old)
    corner_difference=reduced.full[2,2]-kernel.reduced.base[2,2]
    assert corner_difference.real>0
    proof['29']={'source_corner_change_lower':api.b_checks.lower(corner_difference.real)}
    receipt,frames=api.certify_inertia(base,(.25,1.,noise.r))
    assert receipt['status']=='EXCLUDED'
    proof['30']=receipt

    changed=api.MaterialNode('modified','control',box.ell,np.array([.375,.5,.5,.5]))
    trial=dict(t=.25,j=1,weights=weights)
    direct=api.certify_box(changed,old,{},dict(trials=[trial]))
    compressed=api.certify_box(changed,old,{'root':kernel},dict(trials=[dict(trial,kernel='root')]))
    assert direct['status']==compressed['status']=='EXCLUDED'
    assert direct['records'][0]['box_hash']==compressed['records'][0]['box_hash']==changed.box_hash
    assert compressed['records'][0]['eliminated_positive']==2
    proof['31']={'direct':direct,'fixed_kernel':compressed}

    extension=api.extend_traces(kernel,api.local_pair_rows(model,[2],old))
    assert extension['kernel'].reduced.dim==4
    proof['32']={k:v for k,v in extension.items() if k!='kernel'}
    with np.load(api.A/'results/acquisition.npz',allow_pickle=False) as src:
        a=api.a_bounds.Acquisition(src['G'],src['S'],src['B']@src['V'],src['Y'][0,0,0],float(src['radii'][0,0,0]))
        case=src['boxes'][0,0,0,0]
    area=tuple(int(x) for x in api.a_chart.REGION)
    bank=api.prepare_a_chart_kernel(a,0.,1.)
    oracle=api.exact_row_oracle(api.MaterialNode('A-control','A',case[0],case[1]),[area[0]],a,bank,
                                dict(weights=(1.,1.,0.),seconds=.8,max_depth=1,max_leaves=2))
    assert 'proof' in oracle and oracle['proof']['leaves']
    proof['33']={'status':oracle['status'],'leaves':len(oracle['proof']['leaves']),
                 'chart_frame_count':len(oracle['frames'])}

    cover=api.build_task_cover(dict(side=12,shape='single',window=(2,3,2,3),b=.125,
                                    low=.8,cap=2.,mode='core_outer',contract_id='B-control'))
    children=cover['refine'](cover['root'])
    assert len(children)==4 and all(c.parent_id==cover['root'].node_id for c in children)
    assert all(np.all(cover['root'].ell<=c.ell) and np.all(c.u<=cover['root'].u) for c in children)
    proof['34']={'centers':len(cover['root'].candidate_centers),'children':len(children),'proof':cover['coverage']}

    actual=np.full(4,.25)
    y=model.S@np.linalg.solve(np.eye(4)-actual[:,None]*model.G,actual[:,None]*model.Bdict[:,[0]])
    small=api.Acquisition(model,old.V,y,noise,old.receiving_basis_id,old.energies,old.power_gram)
    witness=api.search_witness(box,small,dict(max_evaluations=2))
    assert witness['status']=='WITNESS_FOUND' and witness['witness']['verified_original_square_solve']
    proof['35']={'status':witness['status'],'residual_upper':witness['witness']['residual_upper'],
                 'material_in_box':witness['witness']['material_in_box'],'proposal':witness['search']}

    plan=api.c_model.read_json(api.C/'proposal.json')
    design=next(x for x in plan['source_rows'] if x['scene']=='asymmetric_01' and x['energy_contract']=='added')
    scene,contract=design['scene'],design['energy_contract']
    names=list(design['old_classes'])
    with np.load(api.C/'public_model.npz',allow_pickle=False) as data:
        endpoints=(data[scene+'__lo'],data[scene+'__hi'])
    with np.load(api.C/'old'/f'{scene}__{contract}.npz',allow_pickle=False) as data:
        observations=data['Y']
    active=dict(layout=design['layout'],source_indices=api.c_model.OLD,per_beam_energy=1.,Y=observations,
                noise=api.NoiseContract(api.c_model.high(api.c_model.radius_squared(2).sqrt()),
                                          'archived C radius; disjoint receiver-column measurement errors'))
    envelopes=[api.build_response_envelope(api.MaterialNode(label,'C',l,u),
                     design['old_classes'][label]['chosen_envelope'],active)
               for label,l,u in zip(names,*endpoints,strict=True)]
    assert all(v['status']=='BOUNDED_ENVELOPE' for v in envelopes)
    proof['36']={'bounded_envelopes':len(envelopes),'proofs':[e['proof'] for e in envelopes]}
    action=dict(design['batches']['1']['actions']['conditional_ellipsoid'],layout=design['layout'],
                per_beam_energy=1.)
    sigma=float(api.c_model.COLUMN_RADIUS)
    score=api.score_action(action,envelopes,api.NoiseContract(sigma,'same recorded C absolute error',sigma))
    assert score['status']=='STRICT_EDGE' and score['worst_survivors_upper']==1
    proof['37']={'margin_lower':score['proof']['noisy_separation_gap_lower'],
                 'source_gram_floor':score['source_gram_floor'],'survivors_upper':score['worst_survivors_upper']}

    with tempfile.TemporaryDirectory(prefix='protocol_c_') as home:
        root=Path(home)/'fixed_operator_paper'
        local=root/'experiments/C'
        (root/'results').mkdir(parents=True)
        local.mkdir(parents=True)
        shutil.copyfile(api.PAPER/'results/p0.json',root/'results/p0.json')
        for name in ('model.py','design.py','registration.json','public_model.npz'):
            shutil.copyfile(api.C/name,local/name)
        shutil.copytree(api.C/'old',local/'old')
        fresh=api.select_and_freeze_action(local)
        assert not (local/'acquired').exists()
        for a,b in zip(fresh['proposal']['source_rows'],plan['source_rows'],strict=True):
            for size,block in a['batches'].items():
                assert block['actions']==b['batches'][size]['actions']
        proof['38']={'fresh_proposal_written':True,'future_data_present':False,
                     'frozen_actions_reproduced':True,'scenes':len(fresh['proposal']['source_rows'])}
    recorded=api.select_and_freeze_action(api.C)
    proof['38']['existing_readout_references_actual_proposal']=bool(recorded['proposal_hash'])

    saved=api.append_observation(old,other,old.data_hash)
    assert np.array_equal(saved.Y[:,0],old.Y[:,0]) and np.array_equal(saved.V[:,0],old.V[:,0])
    try:
        api.append_observation(old,other,'wrong_old_identity')
    except ValueError:
        rejected=True
    else:
        rejected=False
    assert rejected
    proof['39']={'old_unchanged':True,'changed_identity_rejected':rejected,'joint_radius_upper':saved.noise.r}
    replay=api.replay_certificate(receipt,dict(model=model,node=box,weights=weights,acquisition=old,frames=frames))
    assert replay['status']=='EXCLUDED' and replay['independent_reassembly']
    proof['40']=replay
    schema=json.loads((Path(__file__).parent/'schema.json').read_text())
    assert set(schema['x-protocol-field-map'])==set(map(str,range(9,26)))
    assert all(row['schema'] in schema['$defs'] and row['sources'] and row['boundary']
               for row in schema['x-protocol-field-map'].values())
    assert set(proof)==set(map(str,range(26,41)))
    return dict(status='PASS',scope='Interface and rejection checks; no new experiment',interval_bits=ctx.prec,
                api_rows=proof,schema_rows=sorted(schema['x-protocol-field-map']),elapsed_seconds=perf_counter()-started),frames


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('selfcheck.json'))
    parser.add_argument('--frames',type=Path,default=Path(__file__).with_name('selfcheck_frames.npz'))
    args=parser.parse_args()
    ctx.prec=128
    ctx.threads=1
    with threadpool_limits(1):
        record,arrays=execute()
    args.output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.frames,**arrays)
    print('PASS: 15 API calls and 17 schema mappings (including protocol line 09).')
