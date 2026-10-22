"""Small true-P1 Frobenius diagnostic; no truth, no new sampling target.

Projector fields, rather than signed eigenvectors, define the vector fields.
Two finite-difference scales and reversed Euler commutators distinguish a
numerical bracket from an assertion that a weak distribution is a manifold.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,Book,Target,BudgetStop,save,write
from p3r2.geometry import material_frame


class ProbeBook(Book):
    def __init__(self,*args,**kwargs):self.pilot_RHS=0;super().__init__(*args,**kwargs)
    def _update(self,event,**fields):
        rhs=int(fields.get('rhs',0))
        if self.pilot_RHS+rhs>2048:raise BudgetStop('Projector-bracket2048RHS cap')
        value=super()._update(event,**fields);self.pilot_RHS+=rhs;return value


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--wait-post-scoring',action='store_true');a=ap.parse_args()
    if a.wait_post_scoring:
        while True:
            p=OUT/'POST_SCORING_CHECKPOINT.json';d=json.loads(p.read_text()) if p.exists() else {}
            if d.get('status')=='STOPPED_ERROR':raise RuntimeError('Previous task needs inspection')
            if d.get('status')=='COMPLETE':break
            time.sleep(15)
    b=ProbeBook('C','C_projector_brackets',max_wall=600);status='COMPLETE';error=None
    t=Target(b,keep_fields=True);results=[]
    try:
        write(b.jobdir/'PROTOCOL.json',dict(parents=['DEV_000','DEV_030'],GT_read=False,
            hypothesis='The rank4 material/prior weak projector need not define an integrable distribution',
            fields='X=P(q)c1/||DgP(q)c1||RMS, Y=P(q)c2/||DgP(q)c2||RMS; c1,c2 fixed from center',
            projector='same rank4 material-normalized generalized frame as TrackC; Euclidean latent projector, gauge independent',
            finite_difference_material_scales=[.001,.002],weak_rank=4,property_coordinates='fixed at paid anchor',
            Gaussian_probe_step_cap=.5,new_RHS_cap=2048,maximum_new_states=48,
            commutator='Euler X,Y,-X,-Y versus Y,X,-Y,-X; half material endpoint difference divided by h² cancels common Euler self-drift',
            interpretation='Differential/finite-step diagnostic at exposed P1 anchors, not posterior sampling or exact-fiber discovery'))
        source=b.jobdir/'source/tools/run_p3r2_projector_brackets.py';source.parent.mkdir(parents=True);source.write_bytes(Path(__file__).read_bytes())
        for parent in ('DEV_000','DEV_030'):
            origin=t.load_cache(OUT/f'A/anchors/{parent}.npz',owner='Historical paid P3L y-only anchor; no GT')
            q0=origin.q.copy();J0=t.jacobian(origin);D0=t.prior.derivative(q0)
            f0=material_frame(J0,D0);seeds=f0['Vw'][:,:2].copy()
            M0=np.r_[D0[:,:64].real,D0[:,:64].imag]/24
            Qw=np.linalg.qr(M0@f0['Vw'],mode='reduced')[0]
            y=t.observation(parent)
            def fields(q,frame):
                D=t.prior.derivative(q)[:,:64];directions=[]
                for seed in seeds.T:
                    v=frame['projector']@seed;norm=np.linalg.norm(D@v)/24
                    if norm<1e-12:raise ValueError('Vanishing material-normalized projector field')
                    directions.append(v/norm)
                return directions
            X,Y=fields(q0,f0)
            save(b.jobdir/parent/'CENTER.npz',q=q0,chi=origin.state.chi,J=J0,D=D0,
                 projector=f0['projector'],weak_material_frame=Qw,X=X,Y=Y,eigenvalues=f0['eigenvalues'])
            def evaluate(q,label):
                if np.linalg.norm(q[:64]-q0[:64])>.5:raise ValueError('Frozen Gaussian probe cap')
                p=t.evaluate(q,parent+'_'+label);J=t.jacobian(p);D=t.prior.derivative(q);frame=material_frame(J,D)
                A,B=fields(q,frame)
                eigenvalues=frame['eigenvalues'];gap=None if len(eigenvalues)<5 else float((eigenvalues[4]-eigenvalues[3])/max(eigenvalues[4],1.))
                save(b.jobdir/parent/(label+'.npz'),q=q,chi=p.state.chi,prediction=p.state.data,
                     J=J,projector=frame['projector'],X=A,Y=B,eigenvalues=eigenvalues,
                     material_rank=frame['material_rank'])
                return dict(q=q.copy(),chi=p.state.chi.copy(),prediction=p.state.data.copy(),X=A,Y=B,
                    material_rank=frame['material_rank'],relative_selection_gap=gap,
                    Phi=t.phi(p,y),endpoint=p.key)
            parent_rows=[]
            for h in (.001,.002):
                if h*max(np.linalg.norm(X),np.linalg.norm(Y))>.1:
                    parent_rows.append(dict(h=h,status='NOT_RUN_PROBE_TOO_LARGE',new_RHS=0));continue
                label='h'+str(h).replace('.','p');start=b.pilot_RHS
                probes={}
                for name,v in (('X',X),('Y',Y)):
                    for sign in (-1,1):
                        q=q0.copy();q[:64]+=sign*h*v
                        probes[(name,sign)]=evaluate(q,label+'_'+name+('plus' if sign>0 else 'minus'))
                bracket=(probes[('X',1)]['Y']-probes[('X',-1)]['Y']-
                         probes[('Y',1)]['X']+probes[('Y',-1)]['X'])/(2*h)
                material_bracket=M0@bracket
                normal=material_bracket-Qw@(Qw.T@material_bracket)
                loops=[]
                for order in ('XY','YX'):
                    q=q0.copy();current=dict(X=X,Y=Y)
                    for k,(field,sign) in enumerate(((order[0],1),(order[1],1),(order[0],-1),(order[1],-1))):
                        q=q.copy();q[:64]+=sign*h*current[field]
                        current=evaluate(q,label+'_'+order+'_leg'+str(k))
                    loops.append(current)
                complex_bracket=D0[:,:64]@bracket
                loop_material=(loops[0]['chi']-loops[1]['chi'])/(2*h*h)
                loop_response=t.white(loops[0]['prediction']-loops[1]['prediction'])/(2*h*h)
                response_prediction=J0[:,:64]@bracket
                residual=float(np.linalg.norm(loop_material-complex_bracket)/max(np.linalg.norm(complex_bracket),1e-30))
                row=dict(h=h,status='COMPLETE',new_RHS=b.pilot_RHS-start,
                    material_bracket_norm=float(np.linalg.norm(material_bracket)),
                    material_normal_bracket_norm=float(np.linalg.norm(normal)),
                    normal_leakage_fraction=float(np.linalg.norm(normal)/max(np.linalg.norm(material_bracket),1e-30)),
                    loop_material_bracket_relative_error=residual,
                    loop_response_bracket_relative_error=float(np.linalg.norm(loop_response-response_prediction)/max(np.linalg.norm(response_prediction),1e-30)),
                    material_tangent_ranks=[p['material_rank'] for p in probes.values()],
                    relative_selection_gaps=[p['relative_selection_gap'] for p in probes.values()],
                    center_material_rank=f0['material_rank'],
                    loop_endpoint_Phi=[p['Phi'] for p in loops],
                    loop_endpoint_material_RMS=[float(np.linalg.norm(p['chi']-origin.state.chi)/24) for p in loops],
                    interpretation='Small-step diagnostic only; no task-scale uncertainty, credible region or exact-fiber certificate')
                save(b.jobdir/parent/(label+'_BRACKET.npz'),bracket=bracket,material_bracket=material_bracket,
                    normal=normal,loop_material_bracket=loop_material,response_prediction=response_prediction,loop_response=loop_response)
                parent_rows.append(row);print(json.dumps(dict(parent=parent,**row)),flush=True)
            if len(parent_rows)==2 and all(r['status']=='COMPLETE' for r in parent_rows):
                with np.load(b.jobdir/parent/'h0p001_BRACKET.npz') as aa:small=aa['material_bracket']
                with np.load(b.jobdir/parent/'h0p002_BRACKET.npz') as aa:large=aa['material_bracket']
                scale_difference=float(np.linalg.norm(small-large)/max(np.linalg.norm(small),1e-30))
            else:scale_difference=None
            results.append(dict(parent=parent,rows=parent_rows,two_scale_material_bracket_relative_difference=scale_difference))
            del origin
        write(b.jobdir/'RESULTS.json',dict(status='COMPLETE',parents=results,new_RHS=b.pilot_RHS,
            new_physics=True,GT_used=False,no_exact_fiber_or_posterior_claim=True))
    except BaseException as exception:status='ERROR';error=repr(exception);raise
    finally:b.close(status,error)


if __name__=='__main__':main()
