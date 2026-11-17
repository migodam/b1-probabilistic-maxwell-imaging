"""Acquisition simulator reads response bank only; no material GT access."""
from common import *
case,arm,label,request,out=sys.argv[1:]
r=np.load(request);a=r['a'];bank=np.load(ROOT/'private/response_bank.npz')
y=bank[case+('_clean' if arm=='clean' else '_shifted')]@a
if arm!='clean':
 seed=int(hashlib.sha256(case.encode()+np.asarray(a,dtype=np.complex128).tobytes()).hexdigest()[:16],16)
 rng=np.random.default_rng(seed);y=y+EPS*np.exp(2j*np.pi*rng.uniform(size=np.shape(y)))
np.savez(out,y=y)
