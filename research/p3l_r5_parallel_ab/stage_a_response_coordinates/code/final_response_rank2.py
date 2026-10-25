"""Local rank2 response cubature in higher-field space; no weak variable in this chart."""
from run_v3_strata import *
F=OUT/'rank2_higher_field';F.mkdir(exist_ok=True)
def main():
    t=Target('A','A1_rank2_final',output_dir=OUT/'physics/A1_rank2_final');model=loadmodel(t,(12,27));rows=[]
    for pi in range(2):
        for ni in range(3):
            name=f'S2_P{pi}_N{ni}'
            with np.load(OUT/f'v2/datasets/{name}.npz') as z:d=z['d'];sig=float(z['sigma'])
            refjson=json.loads((OUT/f'remaining_spaces/{name}_ORIGINAL_13.json').read_text());ref={**refjson,'material_mean':np.array(refjson['material_mean'])[:,0]+1j*np.array(refjson['material_mean'])[:,1],'mean':np.array(refjson['mean']),'cov':np.array(refjson['cov'])}
            x0=ref['mean'];p0=model.full(x0,'rank2_reference_mean_anchor');J0=model.jac(x0,'rank2_actual_two_direction_basis');U,_=np.linalg.qr(J0,mode='reduced');B0=U.T@J0;s0=U.T@(p0['f']-d)/sig;qs={};diagnostics={}
            for n in ([3,5] if pi==0 else [3]):
                zz,ww=quadnodes(n);nodes=[];lw=[];errors=[];fail=0
                for s,w in zip(zz,ww):
                    if t.book.totals()['states']>=32680:raise BudgetStop('final response own soft state stop preserves parent balance')
                    x=x0+np.linalg.solve(B0,sig*(s-s0));ok=False
                    for it in range(8):
                        if np.any(abs(x)>4):break
                        r=model.full(x,'rank2_inverse_full_forward');err=U.T@(r['f']-d)/sig-s
                        if np.linalg.norm(err)<2e-7:ok=True;break
                        x-=np.linalg.solve(B0,sig*err)
                    if not ok:fail+=1;continue
                    J=model.jac(x,'rank2_response_actual_determinant');sg,ld=np.linalg.slogdet(U.T@J)
                    residual=(r['f']-d)/sig;normal=max(0,float(residual@residual-s@s));lp=-.5*x@x-LOGPRIOR_NORM
                    lw.append(np.log(w)+np.log(2*np.pi)+lp-.5*normal+2*np.log(sig)-ld);nodes.append(r);errors.append(float(np.linalg.norm(err)))
                lz=logsumexp(lw);w=np.exp(np.array(lw)-lz);ss=summarize(nodes,w,d,sig,[],lz);qs[str(n)]=ss;diagnostics[str(n)]=dict(inverse_failures=fail,max_inverse_error=max(errors),scope='local rank2 full-response chart, no weak coordinate; global coverage unqualified');report(F/f'{name}_RESPONSE_{n}.json',dict(summary=ss,diagnostics=diagnostics[str(n)]))
            row=dict(observation=name,cross_coordinate_agreement=difference(list(qs.values())[-1],ref),convergence=difference(qs['5'],qs['3']) if '5'in qs else 'NOT_RUN: remaining budget',qualification='NUMERICAL_LOCAL_NOT_CERTIFIED',rank=2,weak_dimension=0,cost=t.book.totals());rows.append(row);report(F/f'{name}_RESULT.json',row);print(json.dumps(plain(row)),flush=True);model.close()
    report(F/'RUN_STATUS.json',dict(status='CLOSED_NUMERICAL_LOCAL',results=rows,cost=t.book.totals()))
if __name__=='__main__':
    try:main()
    except BaseException as e:report(F/'RUN_FAILURE.json',dict(error=repr(e),traceback=traceback.format_exc()));raise
