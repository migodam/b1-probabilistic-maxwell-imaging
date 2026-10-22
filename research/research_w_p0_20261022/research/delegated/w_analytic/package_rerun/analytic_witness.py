"""Research W: exact, nonmagnetic, one-dimensional Maxwell witness.

No inference performance claims. All lengths are in the free-space wavelength
at f0; noise_sigma is the standard deviation of EACH real quadrature.
The e^{-i omega t} convention is used. The finite-layer transfer calculation
includes all internal reflections. Run with Python, NumPy, SciPy, Matplotlib.
"""
from pathlib import Path
import json
import os
import tempfile
import numpy as np
from scipy.stats import norm
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'research_w_mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / 'analytic_results'
OUT.mkdir(exist_ok=True)

def scattering(layers, freq):
    """Return r_left, r_right, t for air/layers/air at physical faces."""
    freq = np.atleast_1d(freq).astype(float)
    M = np.broadcast_to(np.eye(2, dtype=complex), (len(freq), 2, 2)).copy()
    for refractive_index, thickness in layers:
        phase = 2 * np.pi * freq * refractive_index * thickness
        P = np.empty_like(M)
        P[:, 0, 0] = P[:, 1, 1] = np.cos(phase)
        P[:, 0, 1] = 1j * np.sin(phase) / refractive_index
        P[:, 1, 0] = 1j * refractive_index * np.sin(phase)
        M = P @ M
    a, b, c, d = M[:, 0, 0], M[:, 0, 1], M[:, 1, 0], M[:, 1, 1]
    denominator = a + d - b - c
    return (-a-b+c+d)/denominator, (a-b+c-d)/denominator, 2/denominator

def fixed_plane_transmission(layers, freq, plane=1.):
    t = scattering(layers, freq)[2]
    thickness = sum(d for _, d in layers)
    return t * np.exp(2j*np.pi*np.asarray(freq)*(plane-thickness))

def full_target_reflection(layers, freq, position, rho=-1.):
    rl, rr, t = scattering(layers, freq)
    thickness = sum(d for _, d in layers)
    load = rho*np.exp(4j*np.pi*np.asarray(freq)*(position-thickness))
    return rl + t*t*load/(1-rr*load)

def whitened_distance(a, b, sigma):
    # Equivalent to the Euclidean norm of realified data divided by sigma.
    return float(np.linalg.norm(np.asarray(a)-np.asarray(b))/sigma)

def main():
    wall1, wall2 = [(2., .25)], [(3., 1./6.)]
    sigma, nfreq = .01, 9
    h1 = fixed_plane_transmission(wall1, [1.])[0]**2
    h2 = fixed_plane_transmission(wall2, [1.])[0]**2
    r1 = scattering(wall1, [1.])[0][0]
    r2 = scattering(wall2, [1.])[0][0]
    target1, target2 = 1.5, 1.5-1./12.
    y1 = full_target_reflection(wall1, [1.], target1)[0]
    y2 = full_target_reflection(wall2, [1.], target2)[0]
    tests = {
        'max_center_reflection_abs': float(max(abs(r1),abs(r2))),
        'center_task_response_difference_abs': float(abs(h1-h2)),
        'center_task_phase_difference_rad': float(np.angle(h2/h1)),
        'joint_center_data_difference_abs': float(abs(y1-y2)),
        'joint_target_position_difference_wavelength': abs(target1-target2),
    }
    dense = np.linspace(.8, 1.2, 1001)
    energy_errors = []
    for wall in [wall1,wall2,[(1.,.1),(2.3,.12),(1.4,.08)]]:
        rl,rr,t=scattering(wall,dense)
        energy_errors.append(float(np.max(abs(abs(rl)**2+abs(t)**2-1))))
        # Independent check: reciprocal lossless S^* S = I.
        Smat=np.stack([np.stack([rl,t],axis=-1),np.stack([t,rr],axis=-1)],axis=-2)
        energy_errors.append(float(np.max(abs(np.swapaxes(Smat.conj(),-1,-2)@Smat-np.eye(2)))))
    tests['maximum_lossless_unitarity_defect']=max(energy_errors)
    fresnel_errors = []
    for refractive_index, thickness in [wall1[0], wall2[0]]:
        rho = (1-refractive_index)/(1+refractive_index)
        phase = 2*np.pi*dense*refractive_index*thickness
        denominator = 1-rho**2*np.exp(2j*phase)
        reflection = rho*(1-np.exp(2j*phase))/denominator
        transmission = (1-rho**2)*np.exp(1j*phase)/denominator
        r_transfer, _, t_transfer = scattering([(refractive_index, thickness)], dense)
        fresnel_errors.extend([float(np.max(abs(reflection-r_transfer))),
                               float(np.max(abs(transmission-t_transfer)))])
    tests['maximum_independent_fresnel_difference']=max(fresnel_errors)
    assert tests['max_center_reflection_abs'] < 1e-13
    assert abs(tests['center_task_response_difference_abs']-1)<1e-13
    assert abs(tests['center_task_phase_difference_rad']-np.pi/3)<1e-13
    assert tests['joint_center_data_difference_abs']<1e-13
    assert tests['maximum_lossless_unitarity_defect']<1e-12
    assert tests['maximum_independent_fresnel_difference']<1e-12
    rows=[]
    for bandwidth in [0.,.002,.005,.01,.02,.05,.1,.2]:
        frequencies=np.linspace(1-bandwidth/2,1+bandwidth/2,nfreq)
        a=scattering(wall1,frequencies)[0]
        b=scattering(wall2,frequencies)[0]
        d=whitened_distance(a,b,sigma)
        jointd=whitened_distance(full_target_reflection(wall1,frequencies,target1),full_target_reflection(wall2,frequencies,target2),sigma)
        hd=np.linalg.norm(fixed_plane_transmission(wall1,frequencies)**2-fixed_plane_transmission(wall2,frequencies)**2)/np.sqrt(nfreq)
        rows.append({'fractional_bandwidth':bandwidth,'wall_distance':d,'wall_optimal_equal_prior_test_error':float(norm.cdf(-d/2)),'joint_shifted_target_distance':jointd,'task_response_rms_difference':float(hd)})
    payload={'status':'executed_analytic_1d_Maxwell_check_only','noise_sigma_each_real_quadrature':sigma,'number_of_frequencies':nfreq,'wall1_n_d':wall1,'wall2_n_d':wall2,'target1_position':target1,'target2_position':target2,'tests':tests,'bandwidth_sweep':rows,'limitations':['Infinite planar layers, normal incidence, no lateral aperture.','The central equality is single-frequency; bandwidth breaks this particular pair.','Wall-only reflection is an explicit calibration experiment, not assumed available in target-present data.','The joint shifted-target experiment uses only the combined reflection and a PEC plane target.','No 2D or 3D inverse algorithm, posterior sampler, or speedup benchmark was run.']}
    (OUT/'results.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
    import csv
    with (OUT/'bandwidth_sweep.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    xx=np.linspace(0,.4,600)
    axes[0,0].plot(xx,np.where(xx<=.25,4.,1.),color='#156b87',label='W1: n=2, d=0.25')
    axes[0,0].plot(xx,np.where(xx<=1/6,9.,1.),color='#bd603b',label='W2: n=3, d=1/6')
    axes[0,0].set(xlabel='Depth / lambda0',ylabel='Relative permittivity',title='A. Two passive nonmagnetic walls')
    axes[0,0].legend(frameon=False,fontsize=8)
    ff=np.linspace(.95,1.05,501)
    axes[0,1].plot(ff,abs(scattering(wall1,ff)[0]-scattering(wall2,ff)[0]),color='#156b87')
    axes[0,1].axvline(1,color='gray',ls='--',lw=.8)
    axes[0,1].set(xlabel='f / f0',ylabel='Absolute reflection difference',title='B. Exact reflection match only at f0')
    axes[1,0].plot(ff,abs(fixed_plane_transmission(wall1,ff)**2-fixed_plane_transmission(wall2,ff)**2),color='#bd603b')
    axes[1,0].set(xlabel='f / f0',ylabel='Round-trip response difference',title='C. Downstream propagation remains different')
    bw=np.array([x['fractional_bandwidth'] for x in rows[1:]])*100
    axes[1,1].semilogx(bw,[x['wall_distance'] for x in rows[1:]],'o-',color='#156b87',label='Wall-only data')
    axes[1,1].semilogx(bw,[x['joint_shifted_target_distance'] for x in rows[1:]],'s-',color='#bd603b',label='Wall + shifted PEC target')
    axes[1,1].axhline(1,color='gray',ls='--',lw=.8)
    axes[1,1].set(xlabel='Fractional bandwidth (%)',ylabel='Whitened data distance',title='D. Nine frequencies, sigma=0.01 per quadrature')
    axes[1,1].legend(frameon=False,fontsize=8)
    fig.savefig(OUT/'exact_maxwell_witness.pdf')
    fig.savefig(OUT/'exact_maxwell_witness.png',dpi=180)
    plt.close(fig)
    print(json.dumps(payload,indent=2))

if __name__=='__main__':
    main()
