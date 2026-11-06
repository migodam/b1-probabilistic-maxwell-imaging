"""Public downstream development observations from predeclared validation objects."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import json,time,sys
from pathlib import Path
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Model
from research.research_max_struct_r8.src.warm_features import FeatureOperator

def main():
    begin,cpu=time.perf_counter(),time.process_time();out=ROOT/'data_public/warm_development';out.mkdir(parents=True,exist_ok=False)
    rows=json.loads((ROOT/'data_training/validation_manifest.json').read_text())['records'][:8]
    with np.load(ROOT/'data_training/validation.npz') as z:truth=z['m'][:8].copy()
    rng=np.random.default_rng(2026110608);records=[];costs=[]
    for i,row in enumerate(rows):
        cc=row['condition'];opfile=ROOT/'data_public/operators'/f"warm_n24_s{cc['sources']}_a{cc['aperture_fraction']}.npz"
        with np.load(opfile) as z:p={k:z[k] for k in z.files}
        model=Model(p);t,c=time.perf_counter(),time.process_time();clean=model.forward(truth[i].ravel())['y'];sigma=cc['relative_complex_RMS_noise']*np.linalg.norm(clean)/np.sqrt(clean.size)
        Y=clean+sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape));obs=out/(row['case_id']+'.npz')
        np.savez_compressed(obs,Y=Y,sigma=sigma,noise_radius=np.sqrt(chi2.ppf(.99,2*Y.size)/2),grid_shape=[24,24],material_bounds=[.02,20.])
        f=FeatureOperator(p);X,an,fc=f(Y);ff=out/(row['case_id']+'_features.npz');np.savez_compressed(ff,features=X,**an)
        records.append(dict(case_id=row['case_id'],cohort='VALIDATION_DEVELOPMENT',condition=cc,features_file=str(ff.relative_to(REPO)),public_files=[str(opfile.relative_to(REPO)),str(obs.relative_to(REPO))]))
        costs.append(dict(case_id=row['case_id'],LU=1,complex_RHS=cc['sources'],wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,feature_setup=f.setup,features=fc))
    (out/'cases.json').write_text(json.dumps(dict(records=records),indent=2))
    (ROOT/'audit/validation_downstream_generation.json').write_text(json.dumps(dict(status='GENERATED',material_objects=8,new_independent_materials=0,seed=2026110608,records=costs,wall_seconds=time.perf_counter()-begin,cpu_seconds=time.process_time()-cpu,truth_access='EXPLICIT_SUPERVISED_VALIDATION_ONLY'),indent=2))
    print('Prepared 8 validation objects; no test generated.')

if __name__=='__main__':main()
