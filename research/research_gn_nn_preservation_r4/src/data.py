"""Fresh balanced grouped R4 data. Labels and seeds stay evaluator-only."""
from pathlib import Path
import json,time,uuid
from collections import Counter
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from threadpoolctl import threadpool_limits
from research.research_gn_nn_completion_r3.src.data import material_alpha, font_path, read_npz, dump_json
from research.research_z_alg_r4_diversity.src.physics import Model

SEEN=(1,2,3,5,6,8,9)
SHAPES=('primitive_composites','islands','smooth')
FONTS=('Arial Bold.ttf','Georgia Bold.ttf')
OOD_FONTS=('Times New Roman Bold.ttf','Verdana Bold.ttf')

def digit_alpha(rng,digit,font):
    size=192; canvas=Image.new('L',(size,size),0)
    f=ImageFont.truetype(str(font_path(font)),170)
    box=ImageDraw.Draw(canvas).textbbox((0,0),str(digit),font=f)
    glyph=Image.new('L',(box[2]-box[0],box[3]-box[1]))
    ImageDraw.Draw(glyph).text((-box[0],-box[1]),str(digit),font=f,fill=255)
    scale=min(rng.uniform(.57,.78)*size/glyph.height,.70*size/glyph.width)
    glyph=glyph.resize((max(1,round(glyph.width*scale*rng.uniform(.85,1.1))),max(1,round(glyph.height*scale))),Image.Resampling.LANCZOS)
    glyph=glyph.rotate(float(rng.uniform(-24,24)),resample=Image.Resampling.BICUBIC,expand=True)
    canvas.paste(glyph,(int((size-glyph.width)/2+rng.uniform(-.07,.07)*size),int((size-glyph.height)/2+rng.uniform(-.07,.07)*size)))
    return np.asarray(canvas.resize((24,24),Image.Resampling.BOX),float)/255,dict(family='digits',digit=digit,font=font,generator='fresh_transformed_font')

def design():
    rows=[]
    def add(split,group,digits,ndigit,shapes,nshape,ood=False):
        for amp in (.1,.2,.4):
            for digit in digits:
                for repeat in range(ndigit):
                    fonts=OOD_FONTS if ood else FONTS
                    rows.append(dict(split=split,group=group,family='digits',digit=digit,amplitude=amp,font=fonts[(repeat+digit)%2],ood=ood))
            for family in shapes:
                for repeat in range(nshape):
                    rows.append(dict(split=split,group=group,family=family,amplitude=amp,ood=ood))
    add('train','training',SEEN,12,SHAPES,28)
    add('val','validation',SEEN,2,SHAPES,4)
    add('test','new_instance',SEEN,2,SHAPES,2)
    add('test','unseen_class',(4,7),2,('curved_bar',),2)
    add('test','unseen_style',SEEN,1,SHAPES,1,True)
    return rows

def prepare_dataset(output,operators_path,seed=2026110704):
    output=Path(output); started=time.perf_counter()
    if (output/'manifests/evaluation.json').exists(): raise FileExistsError('Preserve existing grouped dataset')
    old=read_npz(operators_path)
    ops={k:old[k] for k in ('G','E0','S','grid_shape','k','centers','sources','receivers','radius') if k in old}
    ops.update(E0=old.get('E0',old.get('U0')),S=old.get('S',old.get('H')),grid_shape=np.array([24,24]))
    if (ops['G'].shape,ops['E0'].shape,ops['S'].shape)!=((576,576),(576,12),(48,576)):raise ValueError('Frozen R4 dimensions required')
    with threadpool_limits(limits=1):
        _,sv,vh=np.linalg.svd(ops['S'],full_matrices=False)
    ops.update(V48=vh.conj().T,S_singular_values=sv)
    output.mkdir(parents=True,exist_ok=True);np.savez_compressed(output/'operators.npz',**ops)
    model=Model(ops);online=[];evaluation=[];parents=[];receipts=[];indices=Counter()
    with threadpool_limits(limits=1):
        for row in design():
            row=row.copy();split=row.pop('split');index=indices[split];indices[split]+=1
            parent_id=f'{split}_p{index:04d}'; stream=[seed,{'train':1,'val':2,'test':3}[split],index]
            rng=np.random.default_rng(np.random.SeedSequence(stream))
            if row['family']=='digits': alpha,meta=digit_alpha(rng,row['digit'],row['font'])
            elif row['family']=='curved_bar':
                canvas=Image.new('L',(192,192),0);d=ImageDraw.Draw(canvas)
                cx,cy=rng.uniform(85,107,2);rx,ry=rng.uniform(48,71,2)
                d.arc((cx-rx,cy-ry,cx+rx,cy+ry),start=int(rng.uniform(10,80)),end=int(rng.uniform(245,345)),fill=255,width=int(rng.uniform(16,29)))
                alpha=np.asarray(canvas.rotate(float(rng.uniform(-25,25)),resample=Image.Resampling.BICUBIC).resize((24,24),Image.Resampling.BOX),float)/255
                meta=dict(family='curved_bar',generator='single_curved_bar')
            else: alpha,meta=material_alpha(rng,row['family'],row['ood'],index)
            amp=row['amplitude'];bg=float(rng.uniform(.08,.15)); m=bg+amp*alpha.reshape(-1)
            membership=dict(seen_image=split=='train',seen_parent=split=='train',seen_class=(row.get('digit') in SEEN if row['family']=='digits' else row['family'] in SHAPES),seen_font=(row.get('font') in FONTS if row['family']=='digits' else True),seen_rule=not(row['ood'] and row['family']!='digits') and row['family']!='curved_bar')
            parent=dict(parent_id=parent_id,split=split,group=row['group'],seed_sequence=stream,background=bg,amplitude=amp,**meta,**membership)
            parents.append(parent);tick=time.perf_counter();state=model.forward(m,label='fresh_balanced_full_vie_generation')
            truth_path=f'truth/{split}/{parent_id}.npz';(output/truth_path).parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(output/truth_path,m=m,j=state['j'],alpha=alpha,amplitude=np.array(amp),amp=np.array(amp))
            receipts.append(dict(parent_id=parent_id,seconds=time.perf_counter()-tick,state_relative_residual=state['state_relative_residual']))
            direction=rng.normal(size=(48,12))+1j*rng.normal(size=(48,12));direction/=np.linalg.norm(direction)
            for noise in (0.,.001):
                obs_id='obs_'+uuid.uuid5(uuid.NAMESPACE_URL,f'gn-nn-r4/{seed}/{parent_id}/{noise}').hex[:20]
                public_path=f'public/{obs_id}.npz';(output/public_path).parent.mkdir(parents=True,exist_ok=True)
                np.savez_compressed(output/public_path,Y=state['y']+noise*np.linalg.norm(state['y'])*direction,noise=np.array(noise))
                public=dict(obs_id=obs_id,parent_id=parent_id,split=split,noise=noise,public_path=public_path)
                online.append(public);evaluation.append(dict(**public,truth_path=truth_path,amplitude=amp,amp=amp,level='low',background=bg,group=row['group'],**meta,**membership))
            if len(parents)%30==0:print(json.dumps(dict(parents=len(parents),seconds=time.perf_counter()-started)),flush=True)
    header=dict(schema='gn-nn-r4-balanced-grouped-v1',operators_path='operators.npz',root='.',parent_grouped_splits=True)
    dump_json(output/'manifests/online.json',dict(**header,records=online))
    dump_json(output/'manifests/evaluation.json',dict(**header,records=evaluation,parents=parents,seed=seed,seen_semantics='membership in TRAIN generator dictionary; image/parent membership in TRAIN; each parent independently generated'))
    receipt=dict(status='GENERATED_FULL_VIE',parents=len(parents),observations=len(online),split_parents=dict(indices),test_groups=dict(Counter(p['group'] for p in parents if p['split']=='test')),wall_seconds=time.perf_counter()-started,counts=dict(model.counts),forward_receipts=receipts,operators_source=str(operators_path),public_keys=['Y','noise'],online_record_keys=['obs_id','parent_id','split','noise','public_path'],truth_keys=['m','j','alpha','amplitude','amp'],gt_exposed_online=False,test_truth_status='SEALED_NOT_SCORED')
    dump_json(output/'GENERATION_RECEIPT.json',receipt)
    return receipt
