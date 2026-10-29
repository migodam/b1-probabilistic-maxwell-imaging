"""Truth-only generator: never imported by a solver. Holdout requires freeze+GO."""
import argparse,json,os,time,uuid
from pathlib import Path
import numpy as np
from scipy.special import hankel1,jv
from scipy.linalg import solve
from PIL import Image,ImageFont,ImageDraw
ROOT=Path(__file__).resolve().parents[1]
FONTS=Path('/System/Library/Fonts/Supplemental')
def dump(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
def geometry():
 c=json.loads((ROOT/'configs/geometry.json').read_text());n=c['grid_shape'][0];h=c['grid_spacing'];a=c['pulse_disk_radius'];k=c['wavenumber']
 q=(np.arange(n)-(n-1)/2)*h;x,y=np.meshgrid(q,q[::-1]);points=np.c_[x.ravel(),y.ravel()]
 circle=lambda count,r,off:np.c_[r*np.cos(np.arange(count)*2*np.pi/count+off),r*np.sin(np.arange(count)*2*np.pi/count+off)]
 sources=circle(c['source_count'],c['source_circle_radius'],c['source_angle_offset']);receivers=circle(c['receiver_count'],c['receiver_circle_radius'],c['receiver_angle_offset'])
 ext=lambda d:.5j*np.pi*k*a*jv(1,k*a)*hankel1(0,k*d)
 d=np.linalg.norm(points[:,None]-points[None,:],axis=2);np.fill_diagonal(d,1.)
 G=ext(d);np.fill_diagonal(G,.5j*np.pi*k*a*hankel1(1,k*a)-1)
 H=ext(np.linalg.norm(receivers[:,None]-points[None,:],axis=2));U0=.25j*hankel1(0,k*np.linalg.norm(points[:,None]-sources[None,:],axis=2))
 return c,dict(G=G,H=H,U0=U0,E0=U0,S=H,centers=points,sources=sources,receivers=receivers,shape=np.array([n,n]),grid_shape=np.array([n,n]),grid_spacing=h,radius=a,k=k,material_bounds=np.array(c['material_bounds']))
def digit(d,font,rng):
 ss=12;side=24*ss;f=ImageFont.truetype(str(FONTS/font),22*ss);bounds=f.getbbox(str(d));w,h=bounds[2]-bounds[0],bounds[3]-bounds[1]
 im=Image.new('L',(w+4*ss,h+4*ss));draw=ImageDraw.Draw(im);stroke=int(rng.uniform(0,.5)*ss)
 draw.text((2*ss-bounds[0],2*ss-bounds[1]),str(d),font=f,fill=255,stroke_width=stroke)
 target_h=int(rng.uniform(18.5,20)*ss);target_w=int(min(18*ss,im.width*target_h/im.height)*rng.uniform(.90,1.10));im=im.resize((target_w,target_h),Image.Resampling.LANCZOS)
 canvas=Image.new('L',(side,side));dx,dy=rng.uniform(-.65,.65,2)*ss;canvas.paste(im,(int((side-im.width)/2+dx),int((side-im.height)/2+dy)))
 angle=float(rng.uniform(-7,7));canvas=canvas.rotate(angle,Image.Resampling.BICUBIC);alpha=np.asarray(canvas.resize((24,24),Image.Resampling.LANCZOS),dtype=float)/255
 bg=float(rng.uniform(.8,1.5));fg=float(rng.uniform(8,18));material=bg+(fg-bg)*alpha
 return material.ravel(),(alpha>=.5).ravel(),dict(font=font,background=bg,foreground=fg,rotation_degrees=angle,shift_pixels=[float(dx/ss),float(dy/ss)],rendered_size=[target_w/ss,target_h/ss],stroke_width_pixels=stroke/ss)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('split',choices=['dev','holdout']);ap.add_argument('--parent-go',action='store_true');args=ap.parse_args()
 if args.split=='holdout' and (not args.parent_go or not (ROOT/'audit/HOLDOUT_FREEZE.json').exists()):raise SystemExit('Holdout blocked: parent GO and audit/HOLDOUT_FREEZE.json required')
 out=ROOT/'public_data'/args.split;truth=ROOT/'private_truth'/args.split
 if out.exists() and list(out.glob('*.npz')):raise SystemExit('Refusing overwrite existing split')
 start=time.perf_counter();c,ops=geometry();setup=time.perf_counter()-start
 if not (ROOT/'public_data/operators.npz').exists():np.savez_compressed(ROOT/'public_data/operators.npz',**ops)
 entropy=int.from_bytes(os.urandom(16),'big');rng=np.random.default_rng(entropy)
 entries=[(2,'Times New Roman Bold.ttf'),(7,'Courier New Bold.ttf')] if args.split=='dev' else [(d,'Arial Rounded Bold.ttf') for d in range(1,10)]
 rng.shuffle(entries);records=[];total=0
 for d,font in entries:
  ident=('dev_' if args.split=='dev' else 'obs_')+uuid.uuid4().hex[:12];m,mask,params=digit(d,font,rng);t=time.perf_counter()
  A=np.eye(len(m),dtype=complex)-m[:,None]*ops['G'];b=m[:,None]*ops['U0'];j=solve(A,b,assume_a='gen',check_finite=False)
  clean=ops['H']@j;elapsed=time.perf_counter()-t;total+=elapsed;noise_level=c['noise_relative_rms'];sigma=float(noise_level*np.sqrt(np.mean(abs(clean)**2)));noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape));Y=clean+noise
  out.mkdir(parents=True,exist_ok=True);truth.mkdir(parents=True,exist_ok=True)
  np.savez_compressed(out/f'{ident}.npz',Y=Y,sigma=sigma,observation_id=np.array(ident))
  np.savez_compressed(truth/f'{ident}.npz',m=m,mask=mask,Y_clean=clean,alpha=(m-params['background'])/(params['foreground']-params['background']))
  records.append(dict(observation_id=ident,digit=d,parameters=params,generation_seconds=elapsed,state_relative_residual=float(np.linalg.norm(A@j-b)/np.linalg.norm(b)),sigma=sigma,clean_rms=float(np.sqrt(np.mean(abs(clean)**2)))))
 dump(truth/'mapping_PRIVATE.json',dict(entropy=str(entropy),records=records))
 dump(ROOT/f'audit/EXPERIMENT_GENERATION_{args.split.upper()}.json',dict(split=args.split,count=len(entries),noise_relative_rms=noise_level,primary_noiseless=args.split=='holdout' and noise_level==0,geometry_setup_seconds=setup,forward_generation_seconds=total,complex_RHS=len(entries)*c['source_count'],LU=len(entries),green_kernel_evaluations=576**2+576*12+48*576,python=os.sys.version,BLAS_thread_env={k:os.getenv(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']},public_ids=[r['observation_id'] for r in records],truth_exposed_to_solver=False))
 print(json.dumps(dict(split=args.split,count=len(entries),setup_seconds=setup,generation_seconds=total,public_ids=[r['observation_id'] for r in records])))
if __name__=='__main__':main()
