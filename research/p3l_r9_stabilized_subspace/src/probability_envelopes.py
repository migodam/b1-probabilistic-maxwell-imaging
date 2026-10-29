from probability_transfer import *
from scipy.optimize import brentq
from functools import lru_cache

@lru_cache(None)
def nodes(n):return leggauss(n)
def cached_gl(n,a,b):
    x,w=nodes(n);return a+(b-a)*(x+1)/2,w*(b-a)/2
# Numerically identical R6 prior integration; cache only standard GL rules.
import probability_transfer as pt
pt.gl=cached_gl

def run():
    if (RES/'envelope_receipt.json').exists():raise RuntimeError('already executed')
    book=Book('envelope');status='FAILED';error=None;allrows={}
    try:
        ph=Physics(book);ph.load_setup()
        for name in ph.obs:
            mode=json.loads((RES/f'{name}_fresh_mode.json').read_text());mode['center']=np.array(mode['center']);mode['C']=np.array(mode['C']);sd=np.sqrt(np.diag(mode['C']));lo=mode['center']-12*sd;hi=mode['center']+12*sd
            audit=json.loads((RES/f'{name}_whole_support_audit.json').read_text());tail_upper=audit['numerical_scaled_omitted_evidence_upper'];reference=json.loads((RES/f'{name}_reduced_GH31.json').read_text());rows=[]
            for level in (8,16,32):
                t=time.perf_counter();cpu=time.process_time();before=dict(book.counts)
                W=[];V=[];ell=[];upper=[];p=[];events=[];cells=[];shift=mode['phi'];phi0=mode['phi']
                with book.span('numerical_task_centered_envelope'):
                    for i in range(level):
                        for j in range(level):
                            x0,x1=lo[0]+(hi[0]-lo[0])*np.array([i,i+1])/level;y0,y1=lo[1]+(hi[1]-lo[1])*np.array([j,j+1])/level
                            a=complex((x0+x1)/2,(y0+y1)/2);maj=ph.majorant_derivative(a,(x1-x0)/2,(y1-y0)/2)
                            assert maj is not None
                            f=ph.reduced(a);rnorm=np.linalg.norm(pack(f-ph.obs[name]))*np.sqrt(2)/SIGMA;e=maj['white_output_radius']
                            ll=float(np.exp(shift-.5*(rnorm+e)**2));uu=float(np.exp(shift-.5*max(0,rnorm-e)**2))
                            rr,wr=cached_gl(3,x0,x1);ii,wi=cached_gl(3,y0,y1);cellweights=[];cellv=[]
                            for r,wr1 in zip(rr,wr):
                                for im,wi1 in zip(ii,wi):
                                    v=complex(r,im);d,_,_=prior(v,96);cellv.append(v);cellweights.append(float(wr1*wi1*d))
                            cellp=sum(cellweights);p.append(cellp);ell.append(ll);upper.append(uu)
                            for v,w in zip(cellv,cellweights):W.append(w);V.append(v);events.append(HOST.real+v.real>=.75)
                            cells.append(dict(parent='retained_rectangle',box=[x0,x1,y0,y1],prior_mass_numerical=cellp,likelihood_lower_scaled_numerical=ll,likelihood_upper_scaled_numerical=uu,white_response_radius=e,prior_moment_rule='GL3 with R6 fibre GL96',physical_reuse_id='probe/transfer_setup.npz',strict_certificate=None))
                    W=np.asarray(W);V=np.asarray(V);ell=np.repeat(np.asarray(ell),9);upper=np.repeat(np.asarray(upper),9);events=np.asarray(events)
                    ZL=float(W@ell);ZU=float(W@upper+tail_upper)
                    tau=HOST.real+V.real
                    def H(c,lower):
                        positive=np.maximum(tau-c,0);negative=np.maximum(c-tau,0)
                        if lower:return float(W@(ell*positive-upper*negative))-tail_upper*max(c-HOST.real,0)
                        return float(W@(upper*positive-ell*negative))+tail_upper*max(HOST.real+RHI-c,0)
                    left=brentq(lambda c:H(c,True),HOST.real,HOST.real+RHI,xtol=1e-14);right=brentq(lambda c:H(c,False),HOST.real,HOST.real+RHI,xtol=1e-14)
                    eventL=float(W@(ell*events));eventU=float(W@(upper*events));compL=float(W@(ell*(~events)));compU=float(W@(upper*(~events)))
                    peL=eventL/(eventL+compU+tail_upper);peU=(eventU+tail_upper)/(eventU+tail_upper+compL)
                    # Naive positive numerator/denominator interval, same cells.
                    naive=[float(W@(ell*tau)/ZU),float((W@(upper*tau)+(HOST.real+RHI)*tail_upper)/ZL)]
                row=dict(case=name,grid_level=level,cells=level**2,likelihood_nodes=level**2,numerical_prior_nodes=len(V),numerical_evidence_enclosure=[ZL,ZU],numerical_omitted_posterior_upper=tail_upper/ZL,numerical_task_centered_mean_interval=[left,right],numerical_naive_mean_interval=naive,numerical_event_interval=[peL,peU],mean_halfwidth=(right-left)/2,event_halfwidth=(peU-peL)/2,mean_stop_old_abs=(right-left)/2<=1e-5,mean_stop_old_SD=(right-left)/2/reference['material_ROI_SD']<=.01,event_stop=(peU-peL)/2<=.002,reference_mean_inside=left<=reference['mean_real_ROI_material']<=right,reference_event_inside=peL<=reference['ROI256_hard_event']<=peU,wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-cpu,cost_counts_delta={k:book.counts[k]-before.get(k,0) for k in book.counts},qualification='NUMERICAL_ENVELOPE_DIAGNOSTIC; Schur/coefficient/rounding and prior-moment error unbounded; these fields are not strict confidence/certificate intervals',strict_evidence_interval=None,strict_mean_interval=None,strict_event_interval=None)
                rows.append(row)
                with (RES/f'{name}_envelope_cells_L{level}.jsonl').open('w') as out:
                    for cell in cells:out.write(json.dumps(plain(cell))+'\n')
                write(RES/f'{name}_task_envelopes.json',rows)
                print(json.dumps(dict(event='ENVELOPE',case=name,level=level,mean_halfwidth=row['mean_halfwidth'],event_halfwidth=row['event_halfwidth'])),flush=True)
            allrows[name]=rows
        write(RES/'task_envelope_results.json',allrows);status='COMPLETE'
    except BaseException as e:error=repr(e);write(AUD/'envelope_failure.json',dict(error=error));raise
    finally:book.close(status,error)
if __name__=='__main__':run()
