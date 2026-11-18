"""Six Maxwell states, two unknown mirror cells, two coherent feeds and three receivers."""
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import root
from scipy.special import j0, j1

from freq_model import FrequencyTM, Fit, ri
from freq_fold import exact_profile
from spectral_moments import PATTERNS, energy_moments

ROOT = Path(__file__).resolve().parents[1]
V = np.array([1., -1.]) / np.sqrt(2.)
STRONG = np.array([[1.], [1.]]) / np.sqrt(2.)
CENTER = np.full(2, 10.5)
SETUP = dict(centers=[[-.18,.04],[.18,.04],[-.29,-.22],[.29,-.22],[-.12,.30],[.12,.30]],
             sources=[[-.55,-.90],[.55,-.90],[-.83,.62],[.83,.62]],
             receivers=[[0.,-1.03],[0.,1.17],[0.,.77]], k=7.4,radius=.047)

class Mirror:
    def __init__(self, turns, fixed_wave=False):
        self.wave=FrequencyTM(SETUP)
        self.turns=np.asarray(turns)
        self.fixed_wave=fixed_wave
        self.gap=np.linalg.norm(np.array(SETUP['sources'])[::2]-np.array(SETUP['sources'])[1::2],axis=1)

    @property
    def factorizations(self):
        return self.wave.factorizations

    def weights(self,delta):
        kd=SETUP['k']*(1+(0. if self.fixed_wave else delta))*self.gap
        phase=2*np.pi*self.turns*delta
        g=j0(kd)
        dphase=2*np.pi*self.turns
        power=2*(1+g*np.cos(phase))
        kp=0. if self.fixed_wave else SETUP['k']
        dp=2*(-j1(kd)*kp*self.gap*np.cos(phase)-g*np.sin(phase)*dphase)
        e=np.exp(1j*phase)
        a=np.stack((np.ones(2),e),axis=1)/np.sqrt(power[:,None])
        ap=np.stack((np.zeros(2),1j*dphase*e),axis=1)/np.sqrt(power[:,None])-a*(dp/(2*power))[:,None]
        return a,ap,power

    def response(self,m,delta=0.,hessian=False,derivative=False):
        mm=np.r_[m,4.5,4.5,6.2,6.2]
        out=self.wave.response(mm,0. if self.fixed_wave else delta,hessian=hessian,derivative=derivative)
        a,ap,_=self.weights(delta)
        def combine(z,weights):
            z=z.reshape((2,2,3)+z.shape[1:])
            val=np.einsum('gs,gsr...->gr...',weights,z)
            return val.reshape((6,)+val.shape[2:])
        y=combine(out[0],a)
        j=combine(out[1][:,:2],a)
        if derivative:
            kp=0. if self.fixed_wave else 1.
            return y,j,kp*combine(out[2],a)+combine(out[0],ap),kp*combine(out[3][:,:2],a)+combine(out[1][:,:2],ap)
        if hessian:
            return y,j,combine(out[2][:,:2,:2],a)
        return y,j

def run():
    tic=perf_counter()
    output={'geometry':SETUP,'known_host':[4.5,4.5,6.2,6.2],'center':CENTER.tolist(),'runs':[]}
    for label,turns in [('equal_feeds',[0,0]),('physical_delays',[1,2]),('fixed_frequency_phase_control',[1,2])]:
        m=Mirror(turns,fixed_wave=label=='fixed_frequency_phase_control')
        probe=CENTER+.6*V
        spread={str(k):float(np.linalg.norm(m.response(probe,k)[0]-m.response(probe[::-1],k)[0]))
                for k in [-.5,-.25,-.1,-.02,0.,.02,.1,.25,.5]}
        y,j,h=m.response(CENTER,hessian=True)
        _,_,_,dj=m.response(CENTER,derivative=True)
        stable=ri(j)@STRONG
        N=np.linalg.svd(stable,full_matrices=True)[0][:,1:]
        T=ri(np.einsum('rij,i,j->r',h,V,V))/2
        normal_c=N.T@T
        derivative=ri(dj@V)
        normal_d=N.T@derivative
        output[label]={'delays_in_reference_periods':turns,'absolute_gaps':spread,'old_singular_values':np.linalg.svd(ri(j),compute_uv=False).tolist(),
                       'normal_curvature_norm':float(np.linalg.norm(normal_c)),'frequency_jet_norm':float(np.linalg.norm(derivative))}
        if label=='equal_feeds':
            continue
        eps=2.e-5
        fd=(m.response(CENTER,eps)[1]-m.response(CENTER,-eps)[1])/(2*eps)
        output[label]['mixed_jet_FD_relative']=float(np.linalg.norm(fd-dj)/np.linalg.norm(dj))
        for gain in ('known','common_prior'):
            C=N.T@np.column_stack((ri(y),ri(1j*y)))
            lam=(np.linalg.norm(C)**2) if gain=='common_prior' else None
            M=np.eye(len(normal_c)) if lam is None else np.linalg.inv(np.eye(len(normal_c))+2*C@C.T/lam)
            A=normal_c@M@normal_c
            D=normal_c@M@normal_d
            E=max(0.,normal_d@M@normal_d-D**2/A)
            den=np.vdot(derivative,derivative).real
            for amplitude in (.2,.6,1.4):
                basis=(CENTER,V,STRONG,amplitude,True)
                for side in (-1,1):
                    truth=CENTER+side*amplitude*V
                    for pattern in ('opposite_edges','two_negative','two_positive'):
                        mu,var,_=energy_moments(pattern)
                        def predicted(Evalue):
                            return amplitude*np.sqrt(A)/(np.sqrt(2*(var*den+mu**2*Evalue))-side*mu*D/(2*np.sqrt(A)))
                        band=predicted(E)
                        x0=-side/2-3*mu*(band/amplitude)*D/(4*A)
                        off,power=PATTERNS[pattern]
                        def residual(z,detail=False):
                            delta=np.asarray(off)*band*z[1]
                            f=Fit(m,delta,[m.response(truth,t)[0] for t in delta],np.sqrt(power),
                                  gain_calibration=(lam,1.) if lam is not None else None)
                            p,g,hh,sg,sh=exact_profile(f,basis,z[0])
                            return (p.tolist(),sg,sh) if detail else np.array([g,hh])/(A*amplitude**4)
                        sol=root(residual,[x0,1.],tol=1.e-8)
                        p,sg,sh=residual(sol.x,True)
                        Fp=np.column_stack([(residual(sol.x+step)-residual(sol.x-step))/.0004 for step in np.eye(2)*.0002])
                        output['runs'].append(dict(experiment=label,pattern=pattern,side=side,a=amplitude,gain=gain,lambda_gain=lam,A=A,D=D,
                             transverse_energy=E,full_derivative_energy=den,predicted_band=band,
                             predicted_omitting_extra_normal_components=predicted(0.),
                             numerical_band=band*sol.x[1],x=sol.x[0],scaled_stationarity=float(np.linalg.norm(residual(sol.x))),
                             strong_eigenvalue=sh,transversality=float(np.linalg.det(Fp)),collision_material=p))
        output[label]['factorizations_including_synthesis']=m.factorizations
    output['wall_s']=perf_counter()-tic
    (ROOT/'results/mirror_spectrum.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    run()
