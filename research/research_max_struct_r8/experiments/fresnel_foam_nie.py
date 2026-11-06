from pathlib import Path
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,argparse
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
import numpy as np
from research.research_max_struct_r8.src.fresnel_adapter import prepare,BoundedMatrixFreeModel
from research.research_max_struct_r8.src.fresnel_nie import run_paper_tsom
from research.research_z_alg_r4_diversity.src.physics import Regularizer
from fresnel_invert import dump

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--attempt',default='');args=ap.parse_args();name='FoamTwinDielTM_8GHz_n60_independent_real_nie'+('_'+args.attempt if args.attempt else '')
 out=ROOT/'RESULTS/fresnel'/name;out.mkdir(exist_ok=False);cpu=time.process_time();wall=time.perf_counter()
 d=prepare(ROOT/'data_local/fresnel/FoamTwinDielTM.exp',frequency=8,side=60,width=.2,matrix_free=True);model=d['model'];model.reg=Regularizer(model.shape,smooth=0.,tv=0.)
 stages=[dict(mf=mf,beta_material=b,paper_delta=delta) for mf,b,delta in [(6,6.,.01),(10,2.,.001),(16,.5,.001)]]
 dump(out/'FROZEN_CONFIG.json',dict(arm='independent_real_nie',frequency_GHz=8,side=60,width_m=.2,source_dependent_DPCS_rank=5,stages=stages,maxiter_perstage=1000,initialization='air0',real_material_box=[0,20],prior='NONE',corrector='NONE',heldout_selection=False,
  qualification='IndependentEq3/Eq14realnonnegativeNIEadaptation; notauthorcode/notqualifiedmatureNIE; realclippingfiniteboxandholdoutacquisitionaredeclareddeparturesfrompaper'))
 try:
  aux=run_paper_tsom(model,d['Y_train'],np.zeros(model.n),stages,rank=5,maxiter=1000);dump(out/'AUXILIARY_RECEIPT.json',aux['receipt']);dump(out/'auxiliary_history.json',aux['history']);np.savez_compressed(out/'paid_incumbent.npz',m=aux['m'])
  full=BoundedMatrixFreeModel(d['public']);full.reg=Regularizer(full.shape,smooth=0,tv=0);st=full.state(aux['m']);pred=st['Y'];rel=lambda mask:float(np.linalg.norm((pred-d['Y'])[mask])/np.linalg.norm(d['Y'][mask]));sh=np.zeros(d['Y'].shape,bool);sh[:,d['held_sources']]=True;sh &=np.isfinite(d['Y']);rh=d['withheld_mask']&~sh
  receipt=dict(status='COMPLETED_INDEPENDENT_REAL_NIE_ADAPTATION',case='FoamTwinDielTM',frequency_GHz=8,side=60,arm='independent_real_nie',initialization='air0',train_relative_complex_residual=rel(d['train_mask']),withheld_all_relative_complex_residual=rel(d['withheld_mask']),withheld_sources_relative_complex_residual=rel(sh),withheld_receivers_relative_complex_residual=rel(rh),training_channels=int(d['train_mask'].sum()),withheld_channels=int(d['withheld_mask'].sum()),state_relative_residual=st['state_relative_residual'],cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,auxiliary=aux['receipt'],physical_cost_training=dict(model.cost),physical_cost_prediction=dict(full.cost),min_contrast=float(aux['m'].min()),max_contrast=float(aux['m'].max()),GT_metrics='UNAVAILABLE',qualification='Separateindependentpaper-equationrealadaptation; nopriorororiginalcorrector; notsamefivearmcontract; notmatureauthorNIE',no_test_refitting=True)
  dump(out/'RECEIPT.json',receipt);np.savez_compressed(out/'derived.npz',m=aux['m'],prediction=pred,train_mask=d['train_mask'],withheld_mask=d['withheld_mask'],centers=d['centers']);print(json.dumps(receipt),flush=True)
 except BaseException as e:dump(out/'FAILED.json',dict(error=repr(e),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost)));raise
if __name__=='__main__':main()
