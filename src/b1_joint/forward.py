"""Unlabelled material adapter to the paid, existing vector-Mie backend."""
from __future__ import annotations
import json
from pathlib import Path
from .material import decode
from ..b1.physics import Forward

class MaterialForward:
    def __init__(self,root,ledger):
        c=json.loads((Path(root)/'configs/protocol.json').read_text());c['nmax']=14
        self.backend=Forward(c,ledger,Path(root)/'research/joint_p2/results/forward_cache')
    def predict(self,x,k=1.6,kind='P2_FULL_VECTOR_MIE',nmax=14,implementation='matching',cache=True):
        try:d=decode(x)
        except Exception as exc:
            def failed_setup():raise exc
            self.backend.ledger.evaluate('SETUP_FAILURE_'+kind,dict(latent=list(map(float,x)),k=float(k),nmax=nmax,kernel_not_started=True),failed_setup)
            raise
        # -1 means no topology label. The radii/eps arguments bypass old H templates.
        return self.backend.run(-1,x[2],x[3],radii=d['radii'],eps_layers=d['eps_layers'],k=float(k),kind=kind,nmax=nmax,implementation=implementation,cache=cache)
