"""Select continuous-band proof inputs using only old-data witnesses and known physics."""
import json
import numpy as np
from scipy.special import j0
from data14 import load_cases, FrequencyTM, cn_radius, ROOT

def main():
    out={};src=ROOT/'results'
    archive=json.loads((src/'pivot_diameter_witnesses.json').read_text())['cases']
    for c in load_cases():
        if c.cohort=='spectral': continue
        m=FrequencyTM(c.geometry)
        pairs=archive[c.case_id]['actions']
        p=np.array(c.geometry['sources']);sp=np.linalg.norm(p[:,None]-p[None,:],axis=-1)
        records=[]
        for name,row in pairs.items():
            w=row['witness'];x,y=np.array(w['materials'])
            za,zb=1+np.array([complex(*v) for v in w['gains_zeta_ri']])
            for frac in (1.,.65,.35):
                for calib in ('original','known_unit'):
                    a=x; b=x+(y-x)*frac
                    g,h=(za,za+(zb-za)*frac) if calib=='original' else (1.,1.)
                    scores=[abs((gain-1)/c.gain_radius)**2+np.linalg.norm(gain*m.response(mat)[0]-c.old)**2/c.additive_radius**2 for gain,mat in ((g,a),(h,b))]
                    if max(scores)>=.9999999: continue
                    fix=[];coh=[];kl=[]
                    for delta in np.linspace(-.02,.02,13):
                        d=g*m.response(a,delta)[0]-h*m.response(b,delta)[0]
                        gram=j0(c.geometry['k']*(1+delta)*sp)
                        n=np.sqrt(max(0.,3*np.vdot(d,np.linalg.solve(gram,d)).real))
                        fix.append(np.linalg.norm(d)/(2*cn_radius(c.sigma,3)))
                        coh.append(n/(2*cn_radius(c.sigma,1)))
                        kl.append(n/c.sigma)
                    records.append(dict(rms=float(np.linalg.norm(a-b)/np.sqrt(6)),fixed=float(max(fix)),coherent=float(max(coh)),gap_sigma=float(max(kl)),calibration=calib,origin=name,fraction=frac,materials=[a.tolist(),b.tolist()],gains_zeta_ri=[[float(np.real(k-1)),float(np.imag(k))] for k in (g,h)],old_scores=scores))
        out[c.case_id]=records
    (src/'band_probe.json').write_text(json.dumps(out,allow_nan=False)+'\n')
    for k,v in out.items():
        good=[r for r in v if r['rms']>1 and r['coherent']<.90]
        if good:
            best=max(good,key=lambda v:v['rms']);print(k,round(best['rms'],2),round(best['coherent'],2),best['calibration'],flush=True)
if __name__=='__main__':main()
