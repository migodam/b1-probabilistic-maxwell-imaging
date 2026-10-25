"""Finite full-Maxwell audit of the C2 plateau seam; no posterior claim.

This is a predeclared mechanism stress path, not a prior-predictive sample.
The physical derivative acts on chi; the latent derivative at the seam is zero.
"""
import os
for key in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from research.p3l_r5_parallel_ab.shared.physics_adapter import Target,pack,save,write
OUT=Path(__file__).resolve().parent/'seam_maxwell'

def main():
    OUT.mkdir(exist_ok=True)
    if (OUT/'RESULTS.json').exists():raise RuntimeError('No implicit paid rerun')
    cfg=dict(role='MECHANISM_STRESS_PATH_NOT_SBC',pair=[0,64],other_coordinates=0,
      property_x64=0.,boundary=-10/13,
      h=[.001,.002,.004,.008,.016,.032,.064,.128],
      fit_indices=[0,1,2,3],fit_target='log raw full-Maxwell distance vs log h',
      scientific_target='separate cubic decoder seam from Maxwell finite remainder',
      truth_used=False,expected_cost=dict(states=9,forward_rhs=72,derivative_rhs=8))
    write(OUT/'FROZEN.json',cfg)
    t=Target('A','parent_seam_probe',output_dir=OUT/'physics')
    x=np.zeros(66);x[0]=cfg['boundary']
    p=t.evaluate(x,purpose='declared_host_seam',retain_reference=True)
    direction=complex(1.3,.105)*t.prior.mask
    rhs=direction[:,None]*p.state.field
    dv=t.inverse_action(p,rhs,purpose='physical_material_JVP_at_seam')
    with t.book.compute('derivative','seam_receiver_action'):
        jf=np.sqrt(2)*pack(t.backend().S@dv)
    f0=np.sqrt(2)*pack(p.prediction);rows=[]
    _,aux=t.decoder(x)
    if np.max(abs(aux['fraction']))>1e-15:raise RuntimeError('Wrong plateau boundary')
    for h in cfg['h']:
        q=x.copy();q[0]+=h
        pp=t.evaluate(q,purpose='declared_seam_h_'+str(h))
        _,a=t.decoder(q);fraction=float(a['fraction'][t.prior.mask][0])
        df=np.sqrt(2)*pack(pp.prediction)-f0
        lin=fraction*jf
        err=float(np.linalg.norm(df-lin));norm=float(np.linalg.norm(df))
        rows.append(dict(h=h,H=fraction,H_over_h_cubed=fraction/h**3,
          material_RMS=float(np.linalg.norm(pp.chi-p.chi)/24),
          raw_response_distance=norm,physical_linear_prediction_distance=float(np.linalg.norm(lin)),
          full_Maxwell_material_remainder=err,relative_physical_remainder=err/norm,
          decoder_coefficient_limit=2.74625,
          backward_residual=float(pp.state.residual)))
        save(OUT/('h_'+str(h)+'.npz'),x=q,chi=pp.chi,prediction=pp.prediction,
             raw_delta=df,material_tangent_prediction=lin)
        del pp
    inds=cfg['fit_indices']
    slope=float(np.polyfit(np.log([rows[i]['h'] for i in inds]),
                 np.log([rows[i]['raw_response_distance'] for i in inds]),1)[0])
    result=dict(status='COMPLETED',config=cfg,rows=rows,loglog_slope=slope,
       physical_derivative_norm=float(np.linalg.norm(jf)),
       latent_seam_derivative_norm=float(np.linalg.norm(t.decoder_jacobian(x)[:,0])),
       level='FULL_DISCRETE_MAXWELL_NUMERICAL_MECHANISM_NOT_PROBABILITY_CERTIFICATE',
       inference='A zero latent Jacobian at the seam coexists with nonzero physical material sensitivity and finite cubic response.',
       limitations=['single predeclared uniform-offset submodel path','not calibrated posterior widths or full66 reachability'])
    write(OUT/'RESULTS.json',result)
    print(json.dumps(dict(status='COMPLETED',slope=slope,physical_Jnorm=result['physical_derivative_norm'])))
if __name__=='__main__':main()
