"""One known-public 50-update wiring run; no GT or full inversion campaign."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
REPO=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r13.src.nfft_control import run_nfft,joint_value_gradient


def main():
    out=Path(__file__).resolve().parent/'nfft_equations_wiring'
    source=REPO/'research/research_max_struct_r13/data_public/development/A__r10_twod_0004.npz'
    wall,cpu=time.perf_counter(),time.process_time()
    with np.load(source,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
    m0=public.pop('m_init');model=MaxwellModel(public);Y=public['Y']
    # Public rank13 was already frozen by the independent public noise-ball
    # rule for this known input in the adapter wiring record, not selected here.
    rank=13
    u,s,vh=svd(model.S,full_matrices=False,check_finite=False)
    JD=vh[:rank].conj().T@((u[:,:rank].conj().T@Y)/s[:rank,None])
    wd=1/np.sum(abs(Y)**2,axis=0);ws=1/np.sum(abs(JD)**2,axis=0)
    reg=lambda m:model.reg.evaluate(m)[:2]
    rng=np.random.default_rng(1301301)
    m=np.clip(m0,.12,19.9);Z=.01*(rng.normal(size=JD.shape)+1j*rng.normal(size=JD.shape))
    dZ=rng.normal(size=JD.shape)+1j*rng.normal(size=JD.shape);dZ/=np.linalg.norm(dZ)
    dm=rng.normal(size=m.shape);dm/=np.linalg.norm(dm)
    value,gZ,gm,parts=joint_value_gradient(model.G,model.E,model.S,Y,JD,Z,m,wd,ws,reg)
    derivatives=[]
    eps=1e-5
    for label,zdir,mdir in [('complex_current',dZ,np.zeros_like(dm)),('real_material',np.zeros_like(dZ),dm)]:
        plus=joint_value_gradient(model.G,model.E,model.S,Y,JD,Z+eps*zdir,m+eps*mdir,wd,ws,reg)[0]
        minus=joint_value_gradient(model.G,model.E,model.S,Y,JD,Z-eps*zdir,m-eps*mdir,wd,ws,reg)[0]
        fd=(plus-minus)/(2*eps);analytic=float(np.vdot(gZ,zdir).real+np.dot(gm,mdir))
        relative=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-12)
        assert relative<1e-5,(label,fd,analytic)
        derivatives.append(dict(direction=label,finite_difference=fd,analytic=analytic,relative_difference=relative))
    audits=[];receipts=[]
    def audit(m,st,value,meta):
        assert abs(value-model.physical_score(m,Y))<1e-14
        audits.append(meta)
    result=run_nfft(model,Y,m0,dict(material_box=[.02,20.],dpcs_rank=rank,
        regularizer=reg,maxiter=50,receipt_path=str(out/'EXECUTION.json'),cost_callback=receipts.append),
        check=lambda:None,audit_callback=audit,audit_every=50)
    assert result['optimizer']['accepted_updates']==50
    assert len(audits)==2
    assert receipts[0]['qualification']=='NOT_ORIGINAL_NFFT_SOM'
    assert np.isfinite(result['final_joint_m']).all()
    np.savez_compressed(out/'WIRING_STATES.npz',m0=m0,selected=result['m'],
        final_joint_m=result['final_joint_m'],final_joint_J=result['final_joint_J'],JD=JD)
    report=dict(status='KNOWN_PUBLIC_50_UPDATE_WIRING_ONLY',qualification='NOT_ORIGINAL_NFFT_SOM',truth_access='NONE',
        source=str(source.relative_to(REPO)),derivative_checks=derivatives,
        derivative_check_extra_cost=dict(receiver_SVD=1,joint_objective_gradient_calls=5),
        optimizer=result['optimizer'],joint_history=result['joint_history'],audits=audits,
        execution_receipt=result['execution_receipt'],validation_CPU_seconds=time.process_time()-cpu,
        validation_wall_seconds=time.perf_counter()-wall,
        scope='One previously exposed public case; max50jointupdates and two new derivative directions; no GT/full campaign')
    (out/'VALIDATION.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],derivatives=derivatives,optimizer=result['optimizer'],
        audits=len(audits),CPU_seconds=report['validation_CPU_seconds']),indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
