"""Analytic controls, not P1 posterior experiments and never new Maxwell RHS."""
import os
os.environ['MPLBACKEND']='Agg'
from pathlib import Path
import json,numpy as np
from scipy.integrate import quad
from scipy.special import gamma,gammainc
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parents[1]
def main():
    sigmas=np.logspace(-1,-12,24);w0=.2;rows=[];maxerr=0.
    for sigma in sigmas:
        integrals={}
        for s in [0,1,3]:
            a=(s+1)/6
            exact=(2*sigma*sigma)**a*gamma(a)*gammainc(a,1/(2*sigma*sigma))/6
            scaled=quad(lambda t:t**s*np.exp(-.5*t**6),0,min(8.,sigma**(-1/3)),epsabs=1e-13,epsrel=1e-12)[0]*sigma**((s+1)/3)
            maxerr=max(maxerr,abs(scaled-exact)/exact)
            integrals[s]=exact
        Z=w0+(1-w0)*integrals[0]
        rows.append(dict(sigma=float(sigma),branch_mass=(1-w0)*integrals[0]/Z,
           regular_task=(1-w0)*integrals[1]/Z,decoder_flat_task=(1-w0)*integrals[3]/Z,
           no_atom_branch_mass=1.,no_atom_regular_task=integrals[1]/integrals[0]))
    slopes={key:float(np.polyfit(np.log(sigmas[-6:]),np.log([r[key] for r in rows][-6:]),1)[0])
      for key in ['branch_mass','regular_task','decoder_flat_task','no_atom_regular_task']}
    (P/'audit/theory_controls').mkdir(parents=True,exist_ok=True)
    result=dict(scope='ANALYTIC_TOY_ONLY_NOT_FULL66; p3 alpha1, prior atom .2 plus Uniform[0,1] branch',
       new_Maxwell_RHS=0,independent_quadrature_relative_error=maxerr,slopes=slopes,rows=rows,
       physical_cubic_regular_material_arc_found=False)
    (P/'audit/theory_controls/RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,ax=plt.subplots(1,2,figsize=(11,4.2))
    for k,label in [('branch_mass','atom + branch: mass'),('regular_task','regular task (s=1)'),('decoder_flat_task','decoder seam task (s=3)')]:
        ax[0].loglog(sigmas,[r[k] for r in rows],label=label)
    ax[0].set(xlabel='Noise scale (toy)',ylabel='Normalized posterior functional');ax[0].legend(fontsize=8)
    t=np.logspace(-4,-1,30)
    ax[1].loglog(t,t**3,label='Decoder seam: material versus parameter')
    ax[1].loglog(t**3,t**3,label='Same seam: response versus material')
    ax[1].loglog(t,t**3,ls=':',lw=3,label='Toy physical cubic: response versus material')
    ax[1].set(xlabel='Parameter / true material displacement as labeled',ylabel='Response / material displacement');ax[1].legend(fontsize=7)
    fig.suptitle('Conditional scaling controls: no Full66 probability claim')
    fig.tight_layout();fig.savefig(P/'figures/THEORY_SCALING_CONTROLS.png',dpi=180);plt.close(fig)
    print(json.dumps(dict(independent_error=maxerr,slopes=slopes,new_RHS=0),indent=2))
if __name__=='__main__':main()
