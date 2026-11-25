"""Independent outward replay and boundary controls, separate from six supplied scripts."""
import json,time
from pathlib import Path
from itertools import product
import numpy as np
from flint import arb,acb,ctx
from fixed_trace import mat,diag,adj,positive_dimension,build_fixed_kernel,update_box,certify_inertia,assemble_border,local_pair_rows
ROOT=Path(__file__).resolve().parents[1]
ctx.prec=128
start=time.perf_counter()
def zero(A):return all(x.contains(0) for x in A.entries())
def save(report):
    (ROOT/'results/p0_core.json').write_text(json.dumps(report,indent=2))
# Interval reconstruction in original coordinates against cell-by-cell constitutive products.
G=mat([[.125+.25j,.25j,0],[.25j,-.25+.25j,.125],[0,.125,.125+.25j]])
S=mat([[1,1j,0],[0,1,-1j]])
B=mat([[1,.5j],[.25,1],[1j,.5]])
V=mat([[1,0,1],[0,1,0]])
X=diag([.25,.5,.75])
J=(diag([1]*3)-X*G).solve(X*B*V)
Y=S*J
T,E=local_pair_rows(G,B,V,range(3))
physical_Z=mat(J.tolist()+diag([1]*3).tolist())
rec=[]
for w,la in (([1,2,1],[0,1,-1]),([0,0,0],[0,0,0]),([1,1,1],[0,0,0])):
    A=assemble_border(G,S,B,V,Y,[0]*3,[1]*3,w,la)
    val=adj(physical_Z)*A*physical_Z
    field=B*V+G*J
    expected=adj(field)*diag([arb(w[i])*x*(x-1) for i,x in enumerate([.25,.5,.75])])*field
    assert zero(val-expected)
    for j in range(1,4):
        receipt=certify_inertia(A,3,arb(1)/256,j,arb(1)/(128*j))
        assert receipt['decision']=='UNRESOLVED'
        rec.append(receipt)
# A congruence leaves the source part fixed, hence its threshold projector unchanged.
transform=mat([[1,0,0,1,0,0],[0,2,0,0,1,0],[0,0,1,0,0,1]]+np.column_stack([np.zeros((3,3)),np.eye(3)]).tolist())
ps=diag([0]*3+[1]*3);assert zero(adj(transform)*ps*transform-ps)
# Full, compressed and exact zero subspaces; near zero is a separate sign test.
a=diag([4,3,2,-1,0,2]);t=mat([[1,1,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
e=mat([[1],[0],[0]]);f=mat([x+y for x,y in zip(t.tolist(),e.tolist())])
ker=build_fixed_kernel(a,t,e,5)
pair=[]
for dyad in ((0,0,0),(1,-2,0),(-2,1,1)):
    du=diag(dyad);full=a+adj(f)*du*f; reduced=update_box(ker,du)
    for threshold in (arb(0),arb(1)/4,arb(2)):
        full_counts=positive_dimension(full-diag([0]*5+[threshold]))
        red_counts=positive_dimension(reduced-diag([0]*ker.dim+[threshold]))
        assert full_counts==[red_counts[0]+ker.eliminated,red_counts[1]]
        pair.append([full_counts,red_counts])
near=positive_dimension(diag([arb(2)**-110,-arb(2)**-111,0]))
assert near==[1,1]
assert certify_inertia(diag([1,1]),1,arb(1),1,arb(1))['decision']=='UNRESOLVED'
# Old measurements with identical errors in multiple columns do not satisfy diag(r_j^2).
assert positive_dimension(mat([[0,-1],[-1,0]]))==[1,1]
# Sharp rank-r relaxation, and source noise at the exact threshold (strict positivity).
h=[6,3,1,-2];tau=[]
for nu in range(5):
    loss=[max(x,0) if i<nu else 0 for i,x in enumerate(h)]
    ds=[max(v-z,0) for v,z in zip(h,loss)]
    assert all(d>=v-z for d,v,z in zip(ds,h,loss)) and sum(ds)==sum(max(x,0) for x in h[nu:])
    tau.append(sum(ds))
# Exhaustive dyadic core sets; nodes split possible centers, not scatterers.
xy=list(product(range(5),repeat=2))
def masks(cs):
    core=[{i for i,x in enumerate(xy) if sum((a-b)**2 for a,b in zip(x,c))<=1} for c in cs]
    outer=[{i for i,x in enumerate(xy) if sum((a-b)**2 for a,b in zip(x,c))<=4} for c in cs]
    return set.intersection(*core),set.union(*outer)
pending=[xy];nodes=0
while pending:
    parent=pending.pop();inner,outer=masks(parent);nodes+=1
    if len(parent)==1:continue
    xs,ys=zip(*parent)
    sx=(min(xs)+max(xs))//2;sy=(min(ys)+max(ys))//2
    children=[[c for c in parent if (c[0]>sx,c[1]>sy)==flag] for flag in product((False,True),repeat=2)]
    assert set().union(*map(set,children))==set(parent)
    for child in filter(None,children):
        a,b=masks(child);assert inner<=a and b<=outer;pending.append(child)
save(dict(status='PASS',precision_bits=128,seconds=time.perf_counter()-start,full_reduced_verified_counts=pair,core_cover_nodes=nodes,near_zero_exact_counts=near,rank_only_sharp_tail=tau,paired_truth_retention=rec))
