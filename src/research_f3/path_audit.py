"""Nine frozen midpoint checks of existing endpoints; no new candidate search."""
import json,traceback,numpy as np
from .common import ROOT,OUT,OLD,Book,dump,rms
from .run import Run,_cpu
from .prior import finite_number

def run():
    cfg=json.loads((OUT/'configs/PATH_VALIDATION.json').read_text());R=object.__new__(Run)
    R.book=Book('F3_PATH_VALIDATION');R.models={};R.pool={};R.jacs={};R.paths={};R.records=[];rows=[];error=None
    try:
        originals={r['case']:r for r in json.loads((OLD/'results/F1/SUMMARY.json').read_text())}
        contexts={}
        for item in cfg['paths']:
            case=item['case']
            if case not in contexts:
                ctx=R.load(originals[case]);f15=json.loads((OUT/f'results/F15/{case}_FINAL.json').read_text());q0=np.array(f15['q'])
                for path in sorted((OUT/'results/endpoints').glob(case+'_F15_*.npz')):
                    with np.load(path) as z:
                        if not np.array_equal(z['q'],q0):continue
                        chi=z['chi'];native=_cpu.FullState(chi,z['current'],z['field'],z['prediction'],np.eye(576)-chi[:,None]*ctx['model'].G,(z['LU'],z['pivots']),float(z['backward_residual']))
                    base=ctx['model']._export_state(native,cached=True,path=path);break
                contexts[case]=(ctx,q0,base)
            ctx,q0,base=contexts[case];end=np.array(item['q'])
            for t in cfg['fractions']:
                q=(1-t)*q0+t*end;key=f'{case}_PATH_{item["id"]}_t{t:g}'
                ep=R.full(ctx,q,key,'PATH_VALIDATION','q_linear_interpolation');phi=R.phi(ctx,ep);lp=ctx['P'].logprior(q)
                rows.append(dict(case=case,path_id=item['id'],parent_method=item['method'],t=t,q=q,array=R.paths[(case,q.tobytes())],
                    material_RMS_from_anchor=rms(ep.chi-base.chi),Phi=phi,delta_Phi=phi-ctx['reference_Phi'],
                    data_feasible=bool(2*phi<=23.785449004057032**2),likelihood_competitive=bool(phi-ctx['reference_Phi']<=2),
                    prior_status='P0_EXACT',log_prior=finite_number(lp),probability_mass=None,
                    role='diagnostic of this fixed parameter path, not optimized continuation or posterior samples'))
                dump(OUT/'results/PATH_RECORDS.json',rows)
    except BaseException as e:error=repr(e);dump(OUT/'PATH_FAILURE.json',dict(error=error,traceback=traceback.format_exc()));raise
    finally:
        R.book.end(error);dump(OUT/'STAGE_CLOSED.json',dict(status='FULL_F3_AUTHORIZATION_CLOSED',error=error,cost=R.book.totals(),no_more_physics=True,pilot_closures_preserved=True))

if __name__=='__main__':run()
