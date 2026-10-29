"""Generate isolated scalar-TM blind targets. Never print materials or seeds."""
import argparse, datetime, importlib.util, json, secrets
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
R2=ROOT.parent/'research_z_alg_r2'
spec=importlib.util.spec_from_file_location('r2_generation_physics',R2/'src/physics.py')
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--split',choices=['dev','holdout','confirm'],required=True)
    ap.add_argument('--count',type=int,default=2);ap.add_argument('--contract',type=Path,required=True)
    args=ap.parse_args(); args.contract=args.contract.resolve();contract=json.loads(args.contract.read_text())
    public=ROOT/'public_data'/args.split; private=ROOT/'private_truth'/args.split
    public.mkdir(parents=True,exist_ok=True);private.mkdir(parents=True,exist_ok=True)
    model=mod.Model(**contract['geometry']);n=model.n
    receipt={'split':args.split,'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'contract':str(args.contract.relative_to(ROOT)),'cases':[],'truth_disclosed':False,'seed_disclosed':False}
    for idx in range(args.count):
        case=f"{contract.get('case_prefix',args.split.upper()+'12')}_{idx:02d}";pub=public/f'{case}.npz';priv=private/f'{case}.npz'
        if pub.exists() or priv.exists():raise FileExistsError(case)
        seed=secrets.randbits(128);rng=np.random.default_rng(seed)
        # Declared prior: two spatial blobs on a known 3x4 disk grid, all values unknown.
        bg=rng.uniform(.7,3.0);m=np.full(n,bg)
        choices=rng.choice(n,2,replace=False)
        for c in choices:
            distance=np.linalg.norm(model.centers-model.centers[c],axis=1)
            support=distance<rng.uniform(.19,.29)
            m[support]=rng.uniform(10.,28.)
        state=model.forward(m,label='blind_generation',independent=True)
        np.savez_compressed(pub,y=state['y'],G=model.G,E0=model.E0,S=model.S,centers=model.centers,sources=model.sources,receivers=model.receivers,radius=model.radius,k=model.k,box=np.asarray(contract['material_box']),case_id=case)
        np.savez_compressed(priv,m=m,seed_decimal=str(seed))
        receipt['cases'].append({'case_id':case,'public_file':str(pub.relative_to(ROOT)),'private_file':str(priv.relative_to(ROOT)),'generation_rhs':model.ns})
    receipt['generation_counts']=dict(model.counts);receipt['generation_timings']=dict(model.timings)
    target=ROOT/'audit'/f'{args.split.upper()}_GENERATION_RECEIPT.json'
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(receipt,indent=2))
    print(json.dumps({'cases':len(receipt['cases']),'generation_rhs':model.ns*args.count,'receipt':str(target),'truth_disclosed':False}))

if __name__=='__main__':main()
