"""One actual measured50-update audit/interrupt wiring check; no heldout or GT."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,time,json,numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.physics_inputs import SourceSpecificPublicModel,run_fresnel_auxiliary
class StopAtPaid50(BaseException):pass
class CheckBudget(BaseException):pass

def main():
    start=time.process_time();path=HERE/'fresnel_public/dielTM_dec8f.npz'
    with np.load(path,allow_pickle=False) as z:public={k:z[k] for k in z.files}
    model=SourceSpecificPublicModel(public);audits=[];callbacks=[]
    def check():
        if time.process_time()-start>10:raise CheckBudget()
    def audit(m,state,value,meta):
        audits.append(dict(**meta,Phi=value,physical_counts=dict(model.cost)))
        if meta.get('audit_kind')=='periodic_material_update':raise StopAtPaid50()
    outcome=None
    try:run_fresnel_auxiliary(model,public['Y'],public['m_init'],'dual_cie',check,audit,
        audit_every=50,cost_callback=callbacks.append)
    except StopAtPaid50:outcome='EXPECTED_STOP_AFTER_PAID_PERIODIC_AUDIT50'
    except CheckBudget:outcome='TEN_CPU_SECOND_CHECK_CAP'
    r=dict(status=outcome,source_public=str(path),truth_access='NONE',heldout_access=False,
        paid_audits=audits,retained_finally_receipt=getattr(model,'r13_auxiliary_cost',None),
        callback_count=len(callbacks),CPU_seconds=time.process_time()-start,
        comparison_or_scientific_result=False)
    (HERE/'FRESNEL_PERIODIC_AUDIT_CHECK.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(dict(status=outcome,CPU=r['CPU_seconds'],audits=len(audits),
        material_updates=r['retained_finally_receipt']['observer_counts']['material_updates_completed'],
        captured_operation_counts=bool(r['retained_finally_receipt']['operation_counts']),callbacks=len(callbacks))))
    if outcome!='EXPECTED_STOP_AFTER_PAID_PERIODIC_AUDIT50':raise RuntimeError('50-update wiring check did not complete within bound')
if __name__=='__main__':main()
