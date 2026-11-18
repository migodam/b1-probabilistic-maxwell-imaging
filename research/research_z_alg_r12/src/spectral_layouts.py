"""Six independently drawn physical layouts, looking for regular interior folds without recovery scores."""
import json
from time import perf_counter

import numpy as np
from scipy.optimize import brentq

from freq_fold import DATA, corrected_pair
from freq_model import FrequencyTM, ri
from spectral_moments import ROOT, investigate


def one(layout_index):
    rng=np.random.default_rng(710430+layout_index)
    original=json.loads(DATA.read_text())
    geometry=dict(original.get('geometry',original))
    xy=np.array(geometry['centers'])
    xy+=rng.uniform(-.055,.055,xy.shape)
    angles=np.array([-2.55,-.53,1.93])+rng.uniform(-.25,.25,3)
    rx=float(rng.uniform(-2.7,2.7)); rad=float(rng.uniform(1.03,1.25))
    geometry.update(centers=xy.tolist(),sources=(1.2*np.stack((np.cos(angles),np.sin(angles)),axis=1)).tolist(),
                    receivers=[[rad*np.cos(rx),rad*np.sin(rx)]],k=float(rng.uniform(5.9,8.8)),
                    radius=float(rng.uniform(.041,.052)))
    model=FrequencyTM(geometry)
    attempts=[]
    # The physical family requires a simple fold. Bracket using sign only, before inversion outcomes.
    for attempt in range(12):
        m0,m1=rng.uniform(3.,25.,(2,len(xy)))
        def signed(t):
            j=ri(model.response(m0+t*(m1-m0))[1])
            return float(np.linalg.det(j/np.linalg.norm(j)))
        l,r=signed(0),signed(1)
        attempts.append({'m0':m0.tolist(),'m1':m1.tolist(),'left':l,'right':r})
        if l*r>=0:continue
        t=brentq(signed,0.,1.,xtol=1e-13)
        center=m0+t*(m1-m0)
        y,j,h=model.response(center,hessian=True)
        u,s,vh=np.linalg.svd(ri(j));v=vh[-1]
        v*=np.sign(v[np.argmax(abs(v))])
        c=float(u[:,-1]@ri(np.einsum('kij,i,j->k',h,v,v)))/2
        if abs(c)<1e-10 or s[-2]<1e-8:continue
        a=.32
        truth=center+a*v
        old=model.response(truth)[0]
        pair=corrected_pair(model,old,center-a*v,center+a*v)
        if np.linalg.norm(pair[1]-pair[0])<.25:
            attempts[-1]['paired_old_solve']='same endpoint'
            continue
        result=investigate(model,pair,f'layout_{layout_index}',(1.,),['two_positive','two_negative','opposite_edges'],True)
        result.update(geometry=geometry,fold_search=attempts,critical_material=center.tolist(),
                      critical_singulars=s.tolist(),curvature=c,nominal_amplitude=a)
        return result
    return {'label':f'layout_{layout_index}','included':False,'geometry':geometry,'fold_search':attempts,
            'reason':'no regular interior sign crossing and retained two-root pair on 12 fixed material chords'}


if __name__=='__main__':
    results=[]
    for idx in range(6):
        result=one(idx);results.append(result)
        (ROOT/'results/spectral_new_layouts.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
        print(result['label'],result['included'],len(result['fold_search']),
              [(x['pattern'],round(x['numerical_band'],5),x['collision_valid'],x['decision']['conditional_TRF_correct_branch']) for x in result.get('runs',[])],flush=True)
