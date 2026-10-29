import json,platform,resource,time
import numpy as np
from scipy.linalg import null_space
from accounting import Book,P,write,identity
from physics import model_for,quality
from charts import MaterialChart,CurrentChart,partition

def run():
    book=Book('INITIAL_PHYSICAL_PROBE');tick=time.perf_counter()
    model=model_for(book);contract=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text())
    m=np.asarray(contract['E1']['public_reference_materials'][0]);state=model.forward(m,jac=True,label='public_reference')
    J=state['J'];dep,free=partition(J);sv=np.linalg.svd(J,compute_uv=False)
    fd=[];rng=np.random.default_rng(2026102901)
    for h in contract['QA']['central_difference_steps']:
        columns=[]
        for i in range(model.n):
            v=np.eye(model.n)[i]
            yp=model.forward(m+h*v,label='FD_plus')['y'];ym=model.forward(m-h*v,label='FD_minus')['y']
            columns.append(model.pack((yp-ym)/(2*h)))
        numerical=np.column_stack(columns)
        fd.append(dict(step=h,relative_error=float(np.linalg.norm(numerical-J)/np.linalg.norm(J))))
    independent=model.forward(m,independent=True,label='independent_numpy')
    with book.op('PERMUTATION_AND_RANK_AUDIT',algebraic_permutations=719):sym=model.permutation_audit()
    material=MaterialChart(model,state['y'],m,dep,free);current=CurrentChart(model,state['y'],state,dep,free)
    z=current.encode(state['j'],m);r,D,_=current.residual_jac(z,m[free])
    zfd=[]
    for h in (1e-4,1e-5,1e-6):
        v=rng.normal(size=len(z));v/=np.linalg.norm(v)
        a=current.residual_jac(z+h*v,m[free],False)[0];b=current.residual_jac(z-h*v,m[free],False)[0]
        zfd.append(dict(step=h,relative_error=float(np.linalg.norm((a-b)/(2*h)-D@v)/np.linalg.norm(D@v))))
    # Receiver-compatible perturbation is NOT a material: test explicit mismatch.
    dark=.03*model.N[:,0,None]*np.array([[1.,1j]])
    jc=state['j']+dark;Ec=model.E0+model.G@jc
    denom=np.sum(abs(Ec)**2,axis=1);hat=np.sum(np.conj(Ec)*jc,axis=1)/denom
    R=jc-hat[:,None]*Ec;minors=jc[:,0]*Ec[:,1]-jc[:,1]*Ec[:,0]
    out=dict(status='PASS',model='asymmetric6',material=m.tolist(),response_singular_values=sv.tolist(),
        numerical_rank={str(c):int(np.sum(sv>c*sv[0])) for c in (1e-8,1e-10,1e-12)},dependent=dep.tolist(),free=free.tolist(),
        source_and_geometry=sym,material_J_FD=fd,current_J_FD=zfd,
        independent_forward=quality(model,independent,state['y'],m),
        current_center_residual=float(np.linalg.norm(r)),current_square_J_singular_values=np.linalg.svd(D,compute_uv=False).tolist(),
        dark_control=dict(data_relative=float(np.linalg.norm(model.S@jc-state['y'])/np.linalg.norm(state['y'])),
            shared_material_residual=float(np.linalg.norm(R)),max_imaginary_recovered_material=float(np.max(abs(hat.imag))),
            normalized_minor_identity_error=float(np.max(abs(np.sum(abs(R)**2,axis=1)-abs(minors)**2/denom))),
            classification='RECEIVER_DARK_BUT_NOT_SHARED_LEGAL_REAL_MATERIAL'),
        first_state_and_J_actual_RHS=3,elapsed_wall_seconds=time.perf_counter()-tick,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    gates=[min(sv)>1e-10*max(sv),max(d['relative_error'] for d in fd)<1e-7,max(d['relative_error'] for d in zfd)<1e-7,
        out['independent_forward']['qualified'],sym['exact_nontrivial_permutations_at_1e_12']==0,np.linalg.norm(R)>1e-5,
        np.linalg.norm(r)<=2e-13,np.linalg.svd(D,compute_uv=False)[-1]>1e-12*np.linalg.svd(D,compute_uv=False)[0]]
    out['status']='PASS' if all(gates) else 'FAIL';out['gates']=list(map(bool,gates))
    write(P/'audit/INITIAL_PHYSICAL_PROBE.json',out)
    np.savez_compressed(P/'audit/PROBE_ARRAYS.npz',m=m,J=J,G=model.G,S=model.S,E0=model.E0,j=state['j'],y=state['y'],N=model.N,dep=dep,free=free,z=z,D=D,dark_j=jc,dark_E=Ec,dark_material=hat,dark_R=R)
    book.end(out['status']);print(json.dumps(out,indent=2))
if __name__=='__main__':run()
