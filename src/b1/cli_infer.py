"""Public observations -> continuous-material event masses and PNG images.
Only frozen paid physics bank/polynomial are reused. No new Maxwell calls.
"""
from __future__ import annotations
import json,time
from pathlib import Path
import numpy as np
from .pilot import dump,compact
from .inference import infer,refinement,top_k,field_summaries,material_image
from .refine import interpolate_bank


def render(output,result,bank,config):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    p=result['topology_probabilities'];order=np.argsort(-np.array(p));fig,axs=plt.subplots(3,2,figsize=(9,11),layout='constrained')
    for row,h in enumerate(order):
        i=result['representative_ids'][h];z,image,void=material_image(h,bank['theta'][i],config)
        for col,values in enumerate([image.real,image.imag]):
            plot=axs[row,col].imshow(values,origin='lower',extent=[z[0],z[-1],z[0],z[-1]],cmap='viridis' if col==0 else 'magma',vmin=1 if col==0 else 0,vmax=2.8 if col==0 else .16)
            axs[row,col].set_title(f'H{h}: {config["topologies"][h]["name"]}\nwhole branch mass {100*p[h]:.4f}% | '+('Re epsilon' if col==0 else 'Im epsilon'))
            fig.colorbar(plot,ax=axs[row,col],shrink=.7)
    fig.suptitle('B1: representative images of continuous topology regions\nExact images have zero mass; out-of-template support is unknown',fontsize=13)
    fig.savefig(output/'top_k_material.png',dpi=160);plt.close(fig)
    fields=field_summaries(result,bank,config);fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,key,title in zip(axs,['void_probability','material_complex_sd'],['Cavity occupancy probability','Complex material standard deviation']):
        im=ax.imshow(fields[key],origin='lower',extent=[-1.1,1.1,-1.1,1.1],cmap='magma');ax.set_title(title);fig.colorbar(im,ax=ax)
    fig.savefig(output/'posterior_fields.png',dpi=160);plt.close(fig)
    np.savez_compressed(output/'material_fields.npz',**{k:v for k,v in fields.items() if isinstance(v,np.ndarray)})


def validate_frozen(root,config):
    frozen=json.loads((Path(root)/'results/DEPLOYMENT_CONTRACT.json').read_text())
    if any(config.get(k)!=v for k,v in frozen.items()):
        raise ValueError('frozen physics/material prior differs from paid bank; cannot reinterpret cached posterior target')


def run(root,config,observation,output,model='refined'):
    root=Path(root);output=Path(output);validate_frozen(root,config)
    if output.exists():raise ValueError('output must be new; keep previous predictions')
    output.mkdir(parents=True);start=time.perf_counter()
    if str(observation).endswith('.npz'):
        a=np.load(observation);y=a['data'];sigma=float(a['noise_real_sd']);pi=a['topology_prior'] if 'topology_prior' in a else np.array(config['topology_prior']);gain_sd=float(a['gain_sd']) if 'gain_sd' in a else 0.
    else:
        a=json.loads(Path(observation).read_text());y=np.array(a['data'],float);sigma=float(a['noise_real_sd']);pi=np.array(a.get('topology_prior',config['topology_prior']));gain_sd=float(a.get('gain_sd',0.))
    if y.shape!=(32,) or not np.all(np.isfinite(y)) or not np.isfinite(sigma) or sigma<=0 or len(pi)!=3 or not np.all(np.isfinite(pi)) or np.any(pi<=0) or abs(pi.sum()-1)>1e-10 or not np.isfinite(gain_sd) or gain_sd<0:raise ValueError('32 finite real measurements, positive fixed noise, normalized proper topology prior required')
    if model=='full-grid':
        bank=dict(np.load(root/'results/banks/order7.npz'));small=dict(np.load(root/'results/banks/order5.npz'));result=infer(bank,y,sigma,pi,gain_sd);q=infer(small,y,sigma,pi,gain_sd);reference=refinement(q,result,config)
        reference['scope']='original full-Mie grid quadrature, no interpolation target'
    elif model=='refined':
        coefficients=np.load(root/'results/refinement/FORWARD_POLYNOMIAL.npz');b17=interpolate_bank(coefficients['coeff7'],17,config,coefficients['bounds']);bank=interpolate_bank(coefficients['coeff7'],25,config,coefficients['bounds']);bcoarse=interpolate_bank(coefficients['coeff5'],17,config,coefficients['bounds'])
        result=infer(bank,y,sigma,pi,gain_sd);r17=infer(b17,y,sigma,pi,gain_sd);c17=infer(bcoarse,y,sigma,pi,gain_sd);q=refinement(r17,result,config);m=refinement(c17,r17,config)
        reference=dict(status='EMPIRICAL_NUMERICAL_PASS' if q['status']==m['status']=='NUMERICAL_REFERENCE_PASS' else 'REFERENCE_UNRESOLVED',integration=q,model_refinement=m,
                       certified_uniform_error=False,independent_validation_evidence='results/refinement/FULL_FORWARD_CHECKS.json; reused development checks, not this observation-specific certificate',
                       target='7x7 polynomial approximation to declared full-Mie model')
    else:raise ValueError('unknown model')
    receipt=dict(model=model,new_maxwell_evaluations=0,measurement_target=dict(physical_dimension='3D_VECTOR',packing='all 16 complex channel real parts then all imaginary parts',real_noise_sd=sigma,complex_noise_sd=float(np.sqrt(2)*sigma),topology_prior=pi,gain_sd=gain_sd),
                 result=compact(result),top_k=top_k(result,bank,config['top_k_display']),reference=reference)
    dump(output/'probabilities.json',receipt);np.savez_compressed(output/'posterior.npz',data=y,noise_real_sd=sigma,topology_prior=pi,gain_sd=gain_sd,posterior_weights=result['weights'],theta=bank['theta'],H=bank['H'],forward=bank['F'])
    render(output,result,bank,config);receipt['wall_seconds_including_images']=time.perf_counter()-start;dump(output/'RUN_RECEIPT.json',receipt)
    print(json.dumps(dict(probabilities=result['topology_probabilities'],remaining_modeled_mass=receipt['top_k']['remaining_modeled_mass'],status=reference['status'],new_physics=0,output=str(output))))
