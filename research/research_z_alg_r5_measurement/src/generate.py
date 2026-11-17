"""Private generator: stdout never contains GT. Runs only after solver freeze."""
from common import *
import secrets
assert (ROOT/'inputs/solver_freeze.json').exists()
book=Ledger('generation');seed=secrets.randbits(128);rng=np.random.default_rng(seed);truth={};bank={};meta={}
anchor=json.loads((R2/'results/general6_GENERAL6_00_general_quotient20.json').read_text())['independently_validated'][0]['m']
for case in ['T0','T1','T2','T3','T4','T5']:
 g=geometry(case);n=len(g['centers']);pos,Q=dictionary(case,g);p=Physics(g,pos,book)
 if case=='T0':m=np.array(anchor)+rng.uniform(-.15,.15,n)
 elif case=='T1':m=rng.uniform(.6,1.5,n)
 elif case=='T5':m=rng.uniform(.5,.7,n)
 else:m=rng.uniform(.5,30,n)
 nominal=p.response(m,np.eye(p.L),independent=True)
 gm=dict(g,radius=.04004);pm=Physics(gm,pos,book);shifted=pm.response(m,np.eye(p.L),independent=True)
 noise=EPS*np.exp(2j*np.pi*rng.uniform(size=p.s));old_clean=nominal[:p.s];old_noisy=shifted[:p.s]+noise
 public=ROOT/'inputs'/f'{case}.npz'
 if public.exists():raise FileExistsError(public)
 np.savez(public,old_clean=old_clean,old_noisy=old_noisy)
 truth[case]=m;bank[case+'_clean']=nominal;bank[case+'_shifted']=shifted
 meta[case]=dict(geometry=g,positions=pos,Q=Q,n=n,k=n,s=p.s,D=6 if n==4 else 20,public_sha=digest(public))
np.savez(ROOT/'private/truth.npz',**truth);np.savez(ROOT/'private/response_bank.npz',**bank)
save(ROOT/'private/seed.json',{'seed':seed});save(ROOT/'inputs/cases.json',meta);save(ROOT/'results/generation_receipt.json',book.receipt())
print('Generated all six targets; GT remains private.')
