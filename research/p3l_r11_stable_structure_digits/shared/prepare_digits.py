from core import *
from PIL import Image,ImageFont,ImageDraw
from matplotlib import font_manager
from scipy.optimize import minimize
from scipy.special import ndtr
import datetime,resource

def main():
 s=load_setup();pr=prior();idx=s['ROI'];start=time.perf_counter()
 if (P/'data/DATA_GENERATED.json').exists():raise RuntimeError('data alreadyfrozen')
 fontpath=font_manager.findfont('DejaVu Sans',fontext='ttf')
 # Font resolved now; every predeclared digit shares same face and pixel size.
 fontpath=str(Path(fontpath).with_name('DejaVuSans-Bold.ttf'))
 font=ImageFont.truetype(fontpath,112)
 cfg=dict(font=fontpath,font_blake2b=blake(fontpath),subpixels=8,physical_grid=24,inkbox_max=[72,104],decoder_fit=dict(objective='mean occupancy squared error +1e-6 norm(z)^2/64',maxiter=1200,ftol=1e-13,gtol=1e-8),representability=dict(occupancy_MSE_max=.025,IoU_min=.75),method='offlinepixeltruthconstructiononly,latentnotavailableonline')
 write(P/'data/GENERATION_CONFIG.json',cfg)
 led=Ledger(P/'data/COST_LEDGER.jsonl',{'phase':'TRUTH_GENERATION'},cap=1024)
 exact=ExactSolver(s['C'],s['O'],s['E'],s['Fh'],ledger=led)
 rows=[];kernel=[];raw=[]
 for d in range(1,10):
  im=Image.new('L',(160,160),0);dr=ImageDraw.Draw(im);bb=dr.textbbox((0,0),str(d),font=font)
  dr.text((-bb[0]+8,-bb[1]+8),str(d),font=font,fill=255)
  box=im.getbbox();glyph=im.crop(box);scale=min(72/glyph.width,104/glyph.height);glyph=glyph.resize((round(glyph.width*scale),round(glyph.height*scale)),Image.Resampling.LANCZOS)
  full=Image.new('L',(192,192),0);full.paste(glyph,((192-glyph.width)//2,(192-glyph.height)//2))
  a=np.asarray(full,dtype=float).reshape(24,8,24,8).mean(axis=(1,3))/255
  a=a.ravel();a[~pr.mask]=0;target=a[idx];raw.append(a)
  def fun(z):
   ss=(pr.config['field_mean']+pr.A[idx]@z)/pr.temperature;aa=saturation(ss);diff=aa-target
   v=np.mean(diff*diff)+1e-6*np.dot(z,z)/64
   g=2*(pr.A[idx].T@(diff*saturation_prime(ss)))/(256*pr.temperature)+2e-6*z/64
   return v,g
  z0=np.linalg.lstsq(pr.A[idx],pr.temperature*(2*target-1),rcond=None)[0]
  opt=minimize(fun,z0,jac=True,method='L-BFGS-B',options=dict(maxiter=1200,ftol=1e-13,gtol=1e-8,maxls=30))
  q=np.r_[opt.x,.5,.5];mat,meta=pr.decode(q);ap=meta['fraction'];diff=ap[idx]-target
  iou=np.sum((ap[idx]>=.5)&(target>=.5))/max(np.sum((ap[idx]>=.5)|(target>=.5)),1)
  status='DECODER_REPRESENTABLE_NUMERIC' if np.mean(diff*diff)<=.025 and iou>=.75 else 'DECODER_REPRESENTABILITY_FAIL'
  for ds,aa,xgt in [('N',a,None),('P',ap,np.r_[opt.x,0.,0.])]:
   oid=f'{ds}_{d-1:03d}';occ=aa[idx];m=np.full(576,HOST);m[idx]+=DELTA_GT*occ
   # Full576 original current equation is truth generator. ROI repeat is independent contractcheck.
   ts=time.perf_counter();M=np.eye(576)-m[:,None]*s['G'];fac=lu_factor(M,check_finite=False);j=lu_solve(fac,m[:,None]*s['Einc'],check_finite=False);F=s['S']@j
   led(dict(event='TRUTH_FULL576',order=576,RHS=8,LU=1,RHS_by_order={'576':8},wall_seconds=time.perf_counter()-ts,parent=oid))
   check=exact.evaluate(occ,DELTA_GT)['F'];wd=np.linalg.norm(white(check-F))
   rng=np.random.default_rng(np.random.SeedSequence([20261030,110011,0 if ds=='N' else 1,d]));eta=SIGMA/np.sqrt(2)*(rng.normal(size=F.shape)+1j*rng.normal(size=F.shape));y=F+eta
   np.savez_compressed(P/'data/online'/f'{oid}.npz',y=y)
   gt=dict(occupancy=aa,material=m,Fmean=F,noise=eta,target_digit=np.array(d),target_occupancy=a,delta=np.array(DELTA_GT))
   if xgt is not None:gt['xgt']=xgt
   np.savez_compressed(P/'data/truth'/f'{oid}.npz',**gt)
   rows.append(dict(parent=oid,dataset=ds,digit=d,noise_seed=[20261030,110011,0 if ds=='N' else 1,d],online_blake2b=blake(P/'data/online'/f'{oid}.npz'),truth_blake2b=blake(P/'data/truth'/f'{oid}.npz'),strict_support=int(np.count_nonzero(occ>0)),truth_noise_phi=phi(F,y),full576_vs_ROI_white=wd,decoder_status=status if ds=='P' else 'PHYSICS_NATIVE_NEW_SUPPORT',projection_MSE=float(np.mean(diff*diff)) if ds=='P' else 0.,projection_IoU=float(iou) if ds=='P' else 1.,latent_norm=float(np.linalg.norm(opt.x)) if ds=='P' else None,fit_success=bool(opt.success),fit_iterations=int(opt.nit)))
   kernel.append(occ)
  print('GENERATED',d,status,'MSE',np.mean(diff*diff),'IoU',iou,'norm',np.linalg.norm(opt.x),flush=True)
 np.savez_compressed(P/'data/kernel_digit_shapes.npz',occupancy=np.stack(kernel),dataset=np.array(['N','P']*9))
 np.savez_compressed(P/'data/digit_targets.npz',occupancy=np.stack(raw))
 write(P/'data/DATA_GENERATED.json',dict(status='COMPLETE',rows=rows,wall_seconds=time.perf_counter()-start,truthgeneration_RHS576=144,ROIcheck_RHS256=144,source_blake2b=blake(Path(__file__)),GT_for_online=False,scope='designed9digits percontract; not randomP1population. P latentsmaybeatypical,diagnosed'))
 led.close()
if __name__=='__main__':main()
