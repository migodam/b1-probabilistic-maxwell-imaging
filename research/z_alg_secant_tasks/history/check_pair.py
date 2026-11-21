"""A single interface check on the two stored roots. No root search."""
import argparse, json
from pathlib import Path
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("historical_results", type=Path)
args=parser.parse_args()
root=args.historical_results
archive=np.load(root/"host_operators.npz")
fibers=json.loads((root/"point_fibers.json").read_text())["fibers"]["minus"]["physical_roots"]
z0,z1=[np.asarray(root["material_midpoint"]) for root in fibers]
zbar,h=(z1+z0)/2,(z1-z0)/2
A=lambda z: np.eye(len(z))-z[:, None]*archive["H"]
current=lambda z: np.linalg.solve(A(z),z[:, None]*archive["B"])
def full(z):
    m=archive["m0"].copy()
    m[archive["active"]]+=z
    return archive["S"]@np.linalg.solve(np.eye(len(m))-m[:,None]*archive["G"],
                                         m[:,None]*archive["E0"])
i=np.linalg.inv(A(zbar))
p=archive["R"]@i
g=archive["H"]@i
e=archive["B"]+archive["H"]@current(zbar)
k=g*h[None,:]
v=np.linalg.solve(np.eye(len(h))-k@k,e)
half=(p*h)@v
drift=(p*h)@k@v
fa,fb,fc=map(full,(z0,z1,zbar))
res={"branch_pair_rms":float(np.sqrt(np.mean((z1-z0)**2))),
     "raw_data_gap": float(np.linalg.norm(fb-fa)),
     "odd_half_identity_error":float(np.linalg.norm(half-(fb-fa)/2)),
     "even_identity_error":float(np.linalg.norm(drift-((fa+fb)/2-fc))),
     "mean_data_vs_mid_material":float(np.linalg.norm(drift)),
     "nominal_raw_upper_noise_radius": float(5*2**-22)}
print(json.dumps(res,indent=2))
