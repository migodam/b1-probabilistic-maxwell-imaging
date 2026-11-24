"""Exposed follow-up: retain all weak source/receiver entry evidence away from one rectangle."""

import argparse
from collections import Counter
import json

import numpy as np
from flint import acb_mat, arb
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from branches import expand
from fast_support import TileProjection
from native_box import HERE,NativeBox,ball,clock,elapsed,frobenius_squared,gershgorin,herm,hi,lo,star
from restricted import score_hard
from site_driver import DEFAULT_PUBLIC
from site_menu import all_sites


def squared_rows(a):
    return [sum((abs(a[i,j])**2 for j in range(a.ncols())),arb(0)) for i in range(a.nrows())]


def positive(x):
    return max(arb(0),x.lower())


class CachedPorts:
    def __init__(self,geometry):
        start=clock()
        self.g=geometry
        self.app=TileProjection(geometry)
        self.app.prepare()
        self.cache={}
        self.cache_spent={}
        for site,branch in all_sites():
            stage=clock()
            ids=np.flatnonzero(expand(branch,geometry.side))
            s,b=geometry.s[:,ids],geometry.b[ids]
            s2,b2=herm(s@s.conj().T),herm(b.conj().T@b)
            _,v=eigh(s2,check_finite=False)
            _,w=eigh(b2,check_finite=False)
            vb,wb=ball(v),ball(w)
            sl,sb=ball(s),ball(b)
            receiver=squared_rows(star(vb)*sl)
            incident=squared_rows(star(sb*wb))
            vnorm=gershgorin(star(vb)*vb,upper=True).sqrt()
            wnorm=gershgorin(star(wb)*wb,upper=True).sqrt()
            entry_caps=[[((a*d).sqrt()*self.app.factor).upper() for d in incident] for a in receiver]
            self.cache[site]=dict(branch=branch,v=vb,w=wb,vw=(vnorm*wnorm).upper(),
                                  caps=entry_caps,receiver=receiver,incident=incident,
                                  receive_gram=sl*star(sl),incident_gram=star(sb)*sb,
                                  receive_gram_numeric=s2,incident_gram_numeric=b2,
                                  scale=float(np.linalg.norm(b)/np.linalg.norm(s)))
            self.cache_spent[site]=elapsed(stage)
        self.preparation=elapsed(start)

    def bound(self,site,y):
        start=clock()
        data=self.cache[site]
        response=star(data['v'])*ball(y)*data['w']
        excess=[]
        for i in range(response.nrows()):
            for j in range(response.ncols()):
                excess.append(positive(abs(response[i,j]).lower()-data['caps'][i][j]))
        value=lo(sum((z*z for z in excess),arb(0)).sqrt()/data['vw'])
        return dict(residual_lower=max(0.,value),significant_entries=sum(z>0 for z in excess),**elapsed(start))

    def block_bound(self,site,y):
        """Distance to the supporting PSD block, verified using negative Ritz values."""
        start=clock()
        data=self.cache[site]
        g=self.g
        a,b=data['receive_gram_numeric'],data['incident_gram_numeric']
        factor=hi(self.app.factor)
        y_ball=ball(y)
        profiles=[]
        best=0.
        for rho in (data['scale']/4,data['scale']/2,data['scale'],2*data['scale'],4*data['scale']):
            hypothesis=herm(np.block([[rho*factor*a,y],[y.conj().T,(factor/rho)*b]]))
            ev,v=eigh(hypothesis,check_finite=False)
            nneg=np.count_nonzero(ev<0)
            if not nneg:
                profiles.append(dict(rho=rho,residual_lower=0.))
                continue
            upper=arb(factor)
            ra=arb(rho)
            ra11=(ra*upper)*data['receive_gram']
            ra22=(upper/ra)*data['incident_gram']
            n=g.nr+g.ns
            matrix=acb_mat(n,n)
            for i in range(n):
                for j in range(n):
                    if i<g.nr and j<g.nr:
                        matrix[i,j]=ra11[i,j]
                    elif i>=g.nr and j>=g.nr:
                        matrix[i,j]=ra22[i-g.nr,j-g.nr]
                    elif i<g.nr:
                        matrix[i,j]=y_ball[i,j-g.nr]
                    else:
                        matrix[i,j]=y_ball[j,i-g.nr].conjugate()
            norm_terms=[]
            for k in range(1,nneg+1):
                w=ball(v[:,:k])
                gap=gershgorin(-star(w)*matrix*w,upper=False)
                if not gap>0:
                    break
                gram=star(w)*w
                if not gershgorin(gram,upper=False)>0:
                    break
                norm_terms.append((gap/gershgorin(gram,upper=True)).lower())
            value=max(0.,lo((sum((z*z for z in norm_terms),arb(0))/2).sqrt()))
            profiles.append(dict(rho=rho,residual_lower=value,negative_eigenvalues_verified=len(norm_terms)))
            best=max(best,value)
        return dict(residual_lower=best,profiles=profiles,**elapsed(start))


def old_fast(row):
    if row['rigor']=='strict_native_translated_tile_port_bound':
        return row['residual_lower']
    return (row.get('translated_tile_check') or {}).get('residual_lower',0.)


def run(args):
    start=clock()
    p=args.public
    online=json.loads(args.online.read_text())
    if online['state']!='SEALED_PUBLIC_SITE_OUTPUTS':
        raise ValueError('This exploratory refinement requires sealed prior site outputs.')
    labelled=json.loads(args.labels.read_text())
    archive=json.loads((p/'manifest.json').read_text())
    models={}
    with np.load(p/'geometry96.npz',allow_pickle=False) as a:
        for ap in ('PW360_10x20','PW180_10x20'):
            g=NativeBox(a['spectrum'],a[ap+'__S'],a[ap+'__B'])
            models[ap]=CachedPorts(g)
    counts=Counter()
    rows={}
    for entry in archive['cases']:
        name=entry['case']
        ref=online['predictions'][name]
        with np.load(p/entry['measurement_file'],allow_pickle=False) as y0:
            y,r=y0['Y'],float(y0['delta'])
        model=models[entry['array']]
        alternatives=[]
        truth=tuple(labelled[entry['object_id']])
        for tile in ref['tile_scores']:
            site=tuple(tile['tile'])
            new=model.bound(site,y)
            old=old_fast(tile['primary'])
            new_prune=new['residual_lower']>r
            direct_prior=tile['primary']['rigor']=='strict_native_exact_support'
            combined=max(old,new['residual_lower'])
            counts['entrywise_prunes']+=new_prune
            counts['prior_rectangle_prunes']+=old>r
            counts['entrywise_avoids_prior_exact_prunes']+=new_prune and direct_prior
            counts['rectangle_prunes_missed_by_entrywise']+=(old>r) and not new_prune
            counts['remaining_after_entrywise_and_rectangle']+=combined<=r
            counts['entrywise_false_prunes']+=site==truth and new_prune
            block=direct=None
            final=combined
            if args.block and combined<=r:
                block=model.block_bound(site,y)
                final=max(final,block['residual_lower'])
                counts['block_extra_prunes']+=final>r
                counts['block_false_prunes']+=site==truth and final>r
            if final<=r and args.replay_fallback:
                branch=model.cache[site]['branch']
                direct=score_hard(model.g,y,r,branch,t_grid=(8.,),alpha_grid=(2.**-14,),max_proofs=1)
                final=max(final,direct['residual_lower'])
                counts['actual_state_fallback_calls']+=1
                counts['actual_state_fallback_prunes']+=direct['decision']=='PRUNE'
            counts['final_false_prunes']+=site==truth and final>r
            alternatives.append(dict(site=site,old_rectangle_lower=old,entrywise=new,block=block,
                                     exact_active=direct,combined_lower=final,decision='PRUNE' if final>r else 'RETAIN'))
        rows[name]=dict(radius=r,tiles=alternatives,stored_truth_label=truth)
        print(name,'entrywise',sum(z['entrywise']['residual_lower']>r for z in alternatives),
              'PSD block extra',sum(z['block'] is not None and z['block']['residual_lower']>r for z in alternatives),flush=True)
    out=dict(scope='Exposed follow-up, only prior sealed compact-site test measurements. No primary decision rule change.',
             geometry_preparation={k:v.preparation for k,v in models.items()},
             port_cache_preparation={k:[dict(tile=site,**cost) for site,cost in v.cache_spent.items()]
                                     for k,v in models.items()},
             counts=dict(counts),observations=rows,**elapsed(start))
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(dict(counts),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--public',type=type(HERE),default=DEFAULT_PUBLIC)
    p.add_argument('--online',type=type(HERE),default=HERE/'site_online.json')
    p.add_argument('--labels',type=type(HERE),required=True)
    p.add_argument('--output',type=type(HERE),default=HERE/'entrywise_exposed.json')
    p.add_argument('--block',action='store_true')
    p.add_argument('--replay-fallback',action='store_true')
    with threadpool_limits(1):
        run(p.parse_args())
