"""Independent narrow qualification of the new toy implementation.

No reconstruction runs. Checks concrete gradient/projection/selection risks.
This script is an audit of software wiring, not a new imaging benchmark.
"""
from pathlib import Path
import importlib.util, json
import numpy as np

P = Path(__file__).parent
spec = importlib.util.spec_from_file_location('toy', P/'consensus_reflection_toy.py')
toy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(toy)

def main():
    xy,G,E,S=toy.setup()
    report={'scope':'Narrow qualification only; no optimization or 3000-step rerun.',
            'checks':{},'cases':[],'total_new_physics_lu':0,'total_new_physics_rhs':0}
    saved=json.loads((P/'toy_results.json').read_text())
    rng=np.random.default_rng(2026110801)
    report['checks']['reciprocity_relative']=float(np.linalg.norm(G-G.T)/np.linalg.norm(G))
    report['checks']['imag_green_min_eigenvalue']=float(np.linalg.eigvalsh((G-G.conj().T)/(2j)).min())
    for row in saved['cases']:
        name=row['case']; a=np.load(P/(name+'.npz'))
        phy=toy.Physics(G,E,S,a['Y'])
        rr={'case':name,'directional_gradient':[],'selection':{},'affine_projection':[]}
        for m in [np.full(len(xy),.1), 2.+.5*rng.random(len(xy))]:
            v=rng.normal(size=len(xy));v/=np.linalg.norm(v)
            f,g=phy.forward(m,True)
            eps=1e-5
            fp=phy.forward(m+eps*v)[0];fm=phy.forward(m-eps*v)[0]
            exact=float(g@v);fd=float((fp-fm)/(2*eps))
            rr['directional_gradient'].append({'analytic':exact,'central_difference':fd,
                'relative_error':abs(exact-fd)/max(abs(exact),abs(fd),1e-12)})
        fzero,gzero=phy.forward(np.zeros(len(xy)),True)
        v=a['truth']/np.linalg.norm(a['truth'])
        fsmall=phy.forward(1e-5*v)[0]
        rr['air_corner']={'f':float(fzero),'gradient_min':float(gzero.min()),
            'gradient_max':float(gzero.max()),'strict_box_KKT_numerical':bool(gzero.min()>1e-10),
            'into_truth_direction_analytic':float(gzero@v),
            'into_truth_direction_one_sided':float((fsmall-fzero)/1e-5)}
        f0,J0,Q0=phy.forward(np.full(len(xy),.1))
        for mode in ['alternating_projection','douglas_rachford']:
            hist=a[mode+'_history'];published=row['methods'][mode]
            fcheck=phy.forward(a[mode])[0]
            best_expected=min(float(f0),float(hist[:,1].min()))
            rr['selection'][mode]={'min_data_objective_in_stored_audits':best_expected,
                'selected_physical_objective':float(fcheck),
                'objective_difference':float(fcheck-best_expected),
                'published_residual_difference':float(np.sqrt(2*fcheck)-published['relative_data_l2'])}
        pr=toy.Projectors(G,E,S,a['Y'],row['delta_CN99'])
        for label,X,Q in [('initial',J0/pr.scale,Q0),
            ('unphysical_random',rng.normal(size=E.shape)+1j*rng.normal(size=E.shape),
             rng.normal(size=E.shape)+1j*rng.normal(size=E.shape))]:
            ax,aq=pr.affine_ball(X,Q)
            residual=pr.T@ax-a['Y']
            gh=(ax-X)+pr.B.conj().T@(aq-Q)
            gd=pr.T.conj().T@residual
            denom=float(np.vdot(gd,gd).real)
            mu=max(0.,float(-np.vdot(gd,gh).real/denom)) if denom>0 else 0.
            stationarity=float(np.linalg.norm(gh+mu*gd)/max(np.linalg.norm(gh),mu*np.linalg.norm(gd),1e-14))
            rr['affine_projection'].append({'input':label,
                'data_bound_ratio':float(np.linalg.norm(residual)/row['delta_CN99']),
                'graph_residual_norm':float(np.linalg.norm(aq-E-pr.B@ax)),
                'estimated_nonnegative_multiplier':mu,
                'relative_KKT_stationarity':stationarity})
        rr['new_physics_lu']=phy.lu_count;rr['new_physics_rhs']=phy.rhs_count
        report['total_new_physics_lu']+=phy.lu_count
        report['total_new_physics_rhs']+=phy.rhs_count
        report['cases'].append(rr)
    # Independent 2x2 symmetric eigensystem reference for joint row projection.
    X=rng.normal(size=E.shape)+1j*rng.normal(size=E.shape)
    Q=rng.normal(size=E.shape)+1j*rng.normal(size=E.shape)
    X[0]=0;Q[1]=0;X[2]=Q[2]=0
    X[3]=1;Q[3]=1j
    X[4]=1e-14;Q[4]=1
    X[5]=1;Q[5]=-1
    cx,cq,m=pr.material(X,Q)
    alpha=m/pr.scale
    actual=np.sum(abs(cx-X)**2+abs(cq-Q)**2,axis=1)
    reference=[]
    for xi,qi in zip(X,Q):
        aa=float(np.vdot(xi,xi).real);bb=float(np.vdot(qi,qi).real)
        cc=float(np.vdot(xi,qi).real)
        vals,vec=np.linalg.eigh(np.array([[aa,-cc],[-cc,bb]]))
        candidates=[0.,pr.amax]
        if abs(vec[0,0])>1e-30:
            alpha_e=float(vec[1,0]/vec[0,0])
            if 0<=alpha_e<=pr.amax:candidates.append(alpha_e)
        reference.append(min((aa-2*cc*z+bb*z*z)/(1+z*z) for z in candidates))
    reference=np.asarray(reference)
    report['checks']['material_projection']={'rows_checked':len(X),
        'max_absolute_objective_difference_to_eigen_reference':float(np.max(abs(actual-reference))),
        'max_relative_objective_difference_to_eigen_reference':float(np.max(abs(actual-reference)/np.maximum(1,abs(reference)))),
        'constitutive_residual_norm':float(np.linalg.norm(cx-alpha[:,None]*cq)),
        'material_box':[float(m.min()),float(m.max())]}
    (P/'toy_qualification_audit.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
