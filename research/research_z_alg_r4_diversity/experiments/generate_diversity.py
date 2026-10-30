"""Local-only full nonlinear VIE generator. Never upload private_truth/."""
import argparse, datetime, json, pathlib, secrets, shutil, sys, time, uuid
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from threadpoolctl import threadpool_limits
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
AMPLITUDES=(.4,2.5,12.)
FONTS={'Arial Bold':'/System/Library/Fonts/Supplemental/Arial Bold.ttf',
       'Georgia Bold':'/System/Library/Fonts/Supplemental/Georgia Bold.ttf'}

def dump(path, value):
    path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def shape_alpha(label):
    """Fixed 384px construction, area average to 24²; no forward calculation."""
    size=384;im=Image.new('L',(size,size));draw=ImageDraw.Draw(im)
    fontname=None
    if label.startswith('digit'):
        digit=int(label[5:]);fontname='Arial Bold' if digit%2 else 'Georgia Bold'
        font=ImageFont.truetype(FONTS[fontname],340)
        box=draw.textbbox((0,0),str(digit),font=font)
        glyph=Image.new('L',(box[2]-box[0],box[3]-box[1]))
        ImageDraw.Draw(glyph).text((-box[0],-box[1]),str(digit),font=font,fill=255)
        # Keep each glyph's aspect ratio, common height and bounded width.
        scale=min(.70*size/glyph.height,.64*size/glyph.width)
        glyph=glyph.resize((round(glyph.width*scale),round(glyph.height*scale)),Image.Resampling.LANCZOS)
        im.paste(glyph,((size-glyph.width)//2,(size-glyph.height)//2))
    elif label=='ring':
        draw.ellipse((65,65,319,319),fill=255);draw.ellipse((123,123,261,261),fill=0)
    elif label=='two_islands':
        draw.ellipse((51,77,163,232),fill=255);draw.ellipse((221,151,333,307),fill=255)
    elif label=='curved_bar':
        draw.arc((53,55,333,329),start=35,end=300,fill=255,width=57)
    else:raise ValueError(label)
    # BOX gives pixel-area alpha; digits/ring preserve antialiased boundaries.
    alpha=np.asarray(im.resize((24,24),Image.Resampling.BOX),float)/255.
    return alpha,fontname

def generate(operators_path):
    from src.physics import Model
    public=ROOT/'public_data';private=ROOT/'private_truth'
    if (public/'manifest.json').exists() or (private/'mapping_PRIVATE.json').exists():
        raise FileExistsError('Preserve existing observations and private mapping')
    with np.load(operators_path,allow_pickle=False) as f:ops={k:f[k] for k in f.files}
    ops['E0']=ops.get('E0',ops.get('U0'));ops['S']=ops.get('S',ops.get('H'));ops['grid_shape']=np.array([24,24]);ops['material_bounds']=np.array([.02,20.])
    if ops['G'].shape!=(576,576) or ops['E0'].shape!=(576,12) or ops['S'].shape!=(48,576):
        raise ValueError('Expected frozen R3 576/12/48 operators')
    if 'k' in ops and not np.isclose(float(ops['k']),2*np.pi):raise ValueError('Expected k=2pi')
    public.mkdir(parents=True,exist_ok=True);private.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(public/'operators.npz',**ops)
    entropy=secrets.randbits(128);rng=np.random.default_rng(entropy)
    model=Model(ops);records=[];generation=[];started=time.perf_counter()
    labels=[f'digit{i}' for i in range(1,10)]+['ring','two_islands','curved_bar']
    with threadpool_limits(limits=1):
        for label in labels:
            alpha,font=shape_alpha(label);background=float(rng.uniform(.08,.15))
            direction=rng.normal(size=(48,12))+1j*rng.normal(size=(48,12))
            direction/=np.linalg.norm(direction)
            for amplitude in AMPLITUDES:
                material=background+amplitude*alpha.ravel()
                before=time.perf_counter();state=model.forward(material,label='generation');clean=state['y']
                born=ops['S']@(material[:,None]*ops['E0'])
                born_full_relative_difference=float(np.linalg.norm(born-clean)/max(np.linalg.norm(clean),1e-300))
                generation.append({'seconds':time.perf_counter()-before,'state_relative_residual':state['state_relative_residual']})
                for noisy in (False,True):
                    ident='obs_'+uuid.uuid4().hex[:16]
                    Y=clean+(1e-3*np.linalg.norm(clean)*direction if noisy else 0.)
                    np.savez_compressed(public/f'{ident}.npz',Y=Y,observation_id=np.array(ident))
                    np.savez_compressed(private/f'{ident}.npz',m=material,alpha=alpha,mask=alpha>=.5,Y_clean=clean)
                    records.append({'observation_id':ident,'shape':label,'kind':'digit' if label.startswith('digit') else 'shape',
                        'font':font,'amplitude':amplitude,'background':background,'noise':'noisy' if noisy else 'clean',
                        'noise_relative':1e-3 if noisy else 0.,'gpu_eligible':label in ('digit2','digit7','ring'),
                        'Born_full_relative_difference':born_full_relative_difference})
    rng.shuffle(records)
    # Publish only opaque IDs and fixed paired-GPU eligibility.
    dump(public/'manifest.json',{'records':[{'observation_id':r['observation_id'],'gpu_eligible':r['gpu_eligible']} for r in records]})
    dump(private/'mapping_PRIVATE.json',{'entropy':str(entropy),'records':records,
        'truth_policy':'LOCAL ONLY; open for scoring only after all predictions sealed',
        'shape_alpha_fixed_across_amplitudes':True,'noise_direction_fixed_per_shape':True})
    receipt={'status':'GENERATED_FULL_VIE','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'observations':len(records),'gpu_observations':sum(r['gpu_eligible'] for r in records),'forward_solves':len(generation),
        'generation_costcounts':dict(model.counts),'generation_wall_seconds':time.perf_counter()-started,
        'forward_receipts':generation,'public_payload_keys':['Y','observation_id'],
        'truth_exposed_to_solver':False,'operator_source':str(operators_path)}
    dump(ROOT/'audit/GENERATION_RECEIPT.json',receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k!='forward_receipts'}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--operators',type=pathlib.Path,default=ROOT.parent/'research_z_alg_r3_digits/public_data/operators.npz')
    ap.add_argument('--go',action='store_true',help='Parent frozen generation GO; writes local private truth')
    args=ap.parse_args()
    if not args.go:ap.error('Generation requires explicit parent frozen GO (--go)')
    generate(args.operators.resolve())
