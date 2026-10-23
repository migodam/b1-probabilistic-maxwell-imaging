"""Paid full-Maxwell probability quadrature bank, with append-only provenance.

Only PUBLIC prior nodes are read. No observations or oracle truth are needed.
CDF numerators integrate the truncated physical interval directly; their
prior interval fraction is kept and never renormalized to one.
"""
import sys,json,gc
from pathlib import Path
import numpy as np
from .accounting import ROOT,Book,write
from .probability_family import Family,CENTERS

class Bank:
    def __init__(self,stage='P1A'):
        self.family=Family.load();self.book=Book(stage);self.root=ROOT/'results/P1A/bank'
        self.root.mkdir(parents=True,exist_ok=True);self.index=self.root/'INDEX.jsonl'
        snapshot=self.root/'PUBLIC_FAMILY_SNAPSHOT.json'
        if snapshot.exists():
            if json.loads(snapshot.read_text())!=self.family.config:raise RuntimeError('cache belongs to a different frozen family/acquisition')
        else:write(snapshot,self.family.config)
        walls=self.root/'PUBLIC_WALLS_SNAPSHOT.npz'
        if walls.exists():
            z=np.load(walls)
            if any(not np.array_equal(z[f'q{i}'],q) for i,q in enumerate(self.family.walls)):raise RuntimeError('cached wall arrays differ from public prior')
        else:np.savez_compressed(walls,q0=self.family.walls[0],q1=self.family.walls[1])
        # Literal source snapshots, not digests: a physics-code change cannot
        # silently reuse an old response as the same computational target.
        files=[ROOT/'src/wr'/n for n in ['physics.py','native_world.py','physical_witness.py','probability_family.py']]
        files += [ROOT.parent/'research_w_continuation/src/wq'/n for n in ['fem.py','meshes.py']]
        files += [ROOT.parent/'research_w_continuation/vendor/p0/geometry.py']
        contract=self.root/'SOURCE_CONTRACT';contract.mkdir(exist_ok=True)
        for i,path in enumerate(files):
            copy=contract/f'{i}_{path.name}'
            if copy.exists():
                if copy.read_bytes()!=path.read_bytes():raise RuntimeError('physics source changed since cache creation: '+str(path))
            else:copy.write_bytes(path.read_bytes())
        self.entries=[json.loads(x) for x in self.index.read_text().splitlines()] if self.index.exists() else []
        self.cache={r['key']:r for r in self.entries}
    def point(self,b,u,c):
        node=self.family.coordinates(b,u,c);b,u,c=node['branch'],node['u'],node['contrast']
        key=f'{int(b)}|{float(u).hex()}|{float(c).hex()}'
        if key in self.cache:
            r=self.cache[key];z=np.load(self.root/r['path']);f=z['F'].copy()
            if (int(z['branch'])!=int(b) or float(z['u'])!=float(u) or float(z['contrast'])!=float(c)
                or f.shape!=(5,4,12) or not np.all(np.isfinite(f))):raise RuntimeError('response cache coordinates/array contract mismatch')
            return f,{**r,'access_incremental_rhs':0,'original_point_rhs':r.get('incremental_rhs',20)}
        number=len(list(self.root.glob('physical_*')))
        dest=self.root/f'physical_{number:05d}'
        data=self.family.point(int(b),float(u),float(c),self.book,str(dest.relative_to(ROOT/'results')))
        r=dict(key=key,branch=int(b),u=float(u),contrast=float(c),path=str((dest/'FULL_RESPONSE.npz').relative_to(self.root)),incremental_rhs=20,frequencies=5)
        with self.index.open('a') as stream:stream.write(json.dumps(r)+'\n');stream.flush()
        self.entries.append(r);self.cache[key]=r
        return data,{**r,'access_incremental_rhs':20,'original_point_rhs':20}

def node_rule(order,cut=None):
    z,w=np.polynomial.legendre.leggauss(order);nodes=[]
    for branch in [0,1]:
        lo=-.01;hi=.01 if cut is None else min(.01,float(cut)-CENTERS[branch][0])
        if hi<=lo:continue
        # density 1/.02 remains the ORIGINAL prior density, including when
        # the CDF interval is truncated. Branch prior remains exactly .5.
        us=(hi+lo)/2+(hi-lo)*z/2;cs=3.+.06*z
        for i,u in enumerate(us):
            for j,c in enumerate(cs):
                mass=.5*(w[i]/2*(hi-lo)/.02)*(w[j]/2)
                nodes.append(dict(branch=branch,u=float(u),contrast=float(c),x=float(CENTERS[branch][0]+u),prior_quadrature_mass=float(mass)))
    return nodes

def main():
    p=json.loads((ROOT/'configs/P1_PROBABILITY_PROTOCOL.json').read_text());bank=Bank()
    rules=[(f'full_{n}',n,None) for n in p['full_orders']]
    rules.extend((f'cdf_{i}_{n}',n,float(cut)) for i,cut in enumerate(p['cdf_checkpoints']) for n in p['reference_orders'])
    out=ROOT/'results/P1A/rules';out.mkdir(parents=True,exist_ok=True)
    for name,order,cut in rules:
        path=out/f'{name}.npz'
        nodes=node_rule(order,cut);data=[];rows=[]
        if path.exists():
            r=np.load(path)
            for name1,name2 in [('branches','branch'),('u','u'),('contrast','contrast'),('weights','prior_quadrature_mass')]:
                if not np.array_equal(r[name1],np.array([n[name2] for n in nodes])):raise RuntimeError('frozen rule changed on resume: '+name)
            if r['F'].shape!=(len(nodes),5,4,12) or not np.all(np.isfinite(r['F'])):raise RuntimeError('invalid saved rule arrays')
            continue
        for node in nodes:
            f,row=bank.point(node['branch'],node['u'],node['contrast']);data.append(f);rows.append({**node,'path':row['path']})
            write(out/'PROGRESS.json',dict(rule=name,completed=len(data),total=len(nodes),paid_unique_points=len(bank.cache)))
            gc.collect()
        np.savez_compressed(path,F=np.array(data),branches=[r['branch'] for r in rows],u=[r['u'] for r in rows],contrast=[r['contrast'] for r in rows],x=[r['x'] for r in rows],weights=[r['prior_quadrature_mass'] for r in rows])
        write(out/f'{name}.json',dict(order=order,cut=cut,nodes=rows,prior_integral=sum(r['prior_quadrature_mass'] for r in rows),probability_object='continuous controlled prior integral approximated by positive quadrature',CDF_measure='direct truncated interval; no bin/node-atom quantiles'))
    write(out/'COMPLETE.json',dict(unique_paid_points=len(bank.cache),complete_rules=[r[0] for r in rules],totals=bank.book.totals()))

if __name__=='__main__':main()
