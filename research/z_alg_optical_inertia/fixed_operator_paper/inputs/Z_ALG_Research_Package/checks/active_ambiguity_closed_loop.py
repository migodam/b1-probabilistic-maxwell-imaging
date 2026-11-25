"""Prospective two-class illumination selection, followed by exact replay.

The score sees only old data and two stated material boxes. It never sees the
new truth response. Selected disk separation and wrong-class exclusion are
replayed with exact rational matrices for exact stored G, B, S.
"""
import json, math
from pathlib import Path
import numpy as np
import sympy as sp

folder=Path(__file__).parent
stored=np.load(folder/'background_compensation_operators.npz')
G,S,Bold=stored['G'],stored['S'],stored['B']
pos=np.array([[0.,0.],[.25,0.],[0.,.25]])
angles=np.arange(16)*2*np.pi/16
Bd=np.exp(1j*2*np.pi*pos@np.array([np.cos(angles),np.sin(angles)]))
old_index=4
Bd[:,old_index]=Bold[:,0]
centers=[np.array([.7,.0625,.0625]),np.array([.74,.03723461520386483,.08442872818623827])]
width=.0005
box_endpoints=[(x-width,x+width) for x in centers]
r_old=r_new=1e-5
t=r_old*r_old
mus=np.logspace(-6,0,7)
# Only the already measured upward datum is supplied to the selection stage.
Yold=S@np.linalg.solve(np.eye(3)-np.diag(centers[0])@G,np.diag(centers[0])@Bold)
Z=np.vstack((np.eye(16)[:,[old_index]],Yold))

def herm(A):return (A+A.conj().T)/2
def predict(lo,hi,mu):
    m=(lo+hi)/2;h=(hi-lo)/2;U=np.eye(3)-np.diag(m)@G
    K=herm(S.conj().T@S+mu*(U.conj().T@U-G.conj().T@np.diag(h*h)@G))
    if np.linalg.eigvalsh(K)[0]<=0:return None
    N=mu*(U.conj().T@np.diag(m)+G.conj().T@np.diag(h*h))@Bd
    L=np.c_[N,S.conj().T]
    D=np.zeros((17,17),complex)
    D[:16,:16]=mu*Bd.conj().T@np.diag(m*m-h*h)@Bd;D[-1,-1]=1
    Q=herm(D-L.conj().T@np.linalg.solve(K,L));H=(Z.conj().T@Q@Z)[0,0].real
    if H>=t:return None
    C=herm(Q-Q@Z@np.linalg.solve(Z.conj().T@Q@Z-t*np.eye(1),Z.conj().T@Q))
    W=C[-1,-1].real
    if W<=0:return None
    cc=-C[-1,:16]/W
    beta=np.abs(C[-1,:16])**2/W-np.diag(C[:16,:16]).real
    if min(beta)<-1e-10:return None
    rr=np.sqrt(np.maximum(0,beta)/W)
    return {'mu':float(mu),'centers':cc,'radii':rr}

options=[]
for lo,hi in box_endpoints:
    options.append([v for mu in mus if (v:=predict(lo,hi,mu)) is not None])
best=(-np.inf,None)
all_best_scores=np.full(16,-np.inf)
for opt0 in options[0]:
    for opt1 in options[1]:
        gaps=np.abs(opt0['centers']-opt1['centers'])-opt0['radii']-opt1['radii']-2*r_new
        gaps[old_index]=-np.inf
        all_best_scores=np.maximum(all_best_scores,gaps)
        j=int(np.argmax(gaps))
        if gaps[j]>best[0]:best=(float(gaps[j]),(j,opt0,opt1))
gap,(selected,opt0,opt1)=best
assert gap>0

# Exact replay starts only after the prospective choice is frozen.
def cq(x):return sp.Rational(float(np.real(x)))+sp.I*sp.Rational(float(np.imag(x)))
Gr=sp.Matrix(3,3,[cq(x) for x in G.ravel()])
Sr=sp.Matrix(1,3,[cq(x) for x in S.ravel()])
Br=sp.Matrix(3,2,[cq(x) for x in Bd[:,[old_index,selected]].ravel()])
xtrue=sp.diag(sp.Rational(7,10),sp.Rational(1,16),sp.Rational(1,16))
Ftrue=sp.simplify(Sr*(sp.eye(3)-xtrue*Gr).inv()*xtrue*Br)
yold=sp.cancel(Ftrue[0,0]);ynew_clean=sp.cancel(Ftrue[0,1])
ro=rn=sp.Rational(1,100000);tr=ro**2
joint_sq=ro**2+rn**2
# Fixed deterministic new error is half its declared radius, along real axis.
ynew=sp.cancel(ynew_clean+rn/2)
Zold=sp.Matrix([1,0,yold])

def hsym(A):return sp.simplify((A+A.conjugate().T)/2)
def exact_predict(box,mu_float):
    lo,hi=[sp.Matrix([sp.Rational(float(v)) for v in end]) for end in box]
    m=(lo+hi)/2;h=(hi-lo)/2
    M=sp.diag(*m);A=sp.diag(*[v*v for v in h]);U=sp.eye(3)-M*Gr
    # Selected weight is a decimal power; preserve its intended exact value.
    mu=sp.Rational(str(float(mu_float)))
    K=hsym(Sr.conjugate().T*Sr+mu*(U.conjugate().T*U-Gr.conjugate().T*A*Gr))
    minors=[sp.simplify(K[:i,:i].det()) for i in range(1,4)]
    assert all(v>0 for v in minors)
    N=mu*(U.conjugate().T*M+Gr.conjugate().T*A)*Br
    L=N.row_join(Sr.conjugate().T)
    D=sp.zeros(3);D[:2,:2]=mu*Br.conjugate().T*(M*M-A)*Br;D[2,2]=1
    Q=hsym(D-L.conjugate().T*K.inv()*L)
    H=sp.simplify((Zold.conjugate().T*Q*Zold)[0]);assert H<tr
    C=hsym(Q-Q*Zold*(H-tr)**-1*Zold.conjugate().T*Q)
    W=sp.simplify(C[2,2]);assert W>0
    cc=sp.cancel(-C[2,1]/W)
    beta=sp.simplify(C[1,2]*C[2,1]/W-C[1,1]);assert beta>=0
    rad_sq=sp.cancel(beta/W)
    Zaug=sp.Matrix([[1,0],[0,1],[yold,ynew]])
    Haug=hsym(Zaug.conjugate().T*Q*Zaug)
    trace=sp.simplify(sp.trace(Haug))
    return {'center':cc,'rad_sq':rad_sq,'trace':trace,'positive_K_minors':minors,'H_old':H,'W':W}

exact=[exact_predict(box_endpoints[0],opt0['mu']),exact_predict(box_endpoints[1],opt1['mu'])]

def sqrt_bounds(q,digits=28):
    q=sp.cancel(q);assert q>=0
    num,den=map(int,sp.fraction(q));scale=10**digits
    floor=math.isqrt(num*scale*scale//den)
    lower=sp.Rational(floor,scale);upper=sp.Rational(floor+1,scale)
    assert lower*lower<=q<=upper*upper
    return lower,upper

dc=sp.expand_complex(exact[1]['center']-exact[0]['center'])
dist_sq=sp.simplify(sp.re(dc)**2+sp.im(dc)**2)
dl,du=sqrt_bounds(dist_sq)
rr=[sqrt_bounds(e['rad_sq']) for e in exact]
gap_lower=dl-rr[0][1]-rr[1][1]-2*rn
assert gap_lower>0
wrong_margin=sp.cancel(exact[1]['trace']-joint_sq)
assert wrong_margin>0
assert (rn/2)**2<joint_sq
baseline_inside=all(sp.Rational(float(lo))<=xtrue[i,i]<=sp.Rational(float(hi))
                    for i,(lo,hi) in enumerate(zip(*box_endpoints[0])))
assert baseline_inside
# Prior exact alternative root remains inside the second interval box.
prior=json.loads((folder/'background_compensation_exact.json').read_text())['proofs'][0]
alt_intervals=[[sp.Rational(37,50),sp.Rational(37,50)]]
alt_intervals += [[sp.Rational(v) for v in prior[key]] for key in ['background1_interval','background2_interval']]
alternative_inside=all(sp.Rational(float(box_endpoints[1][0][i]))<=a<=b<=sp.Rational(float(box_endpoints[1][1][i]))
                       for i,(a,b) in enumerate(alt_intervals))
assert alternative_inside

result={
 'status':'EXACT_STORED_OPERATOR_PROSPECTIVE_SEPARATION_AND_CLOSED_LOOP',
 'scope':'Mechanism example; 3 cells, one receiver; exact rational replay on stored binary64 operators. PSD nu=0, hence not evidence of indefinite advantage or benchmark superiority.',
 'old_source_index':old_index,'old_source_degrees':90,
 'candidate_count':15,'candidate_dictionary_total':16,
 'selected_source_index':selected,'selected_source_degrees':float(angles[selected]*180/np.pi),
 'class_centers':[v.tolist() for v in centers], 'class_halfwidth':width,
 'old_radius':float(ro),'new_radius':float(rn),'joint_radius':float(sp.sqrt(joint_sq)),
 'selected_sector_weights':[opt0['mu'],opt1['mu']],
 'floating_candidate_scores':[None if not np.isfinite(v) else float(v) for v in all_best_scores],
 'certified_clean_disk_centers':[[float(sp.re(e['center'])),float(sp.im(e['center']))] for e in exact],
 'certified_disk_radius_upper':[float(v[1]) for v in rr],
 'certified_noisy_separation_lower':float(gap_lower),
 'certified_noisy_separation_lower_rational':str(gap_lower),
 'new_true_clean_response':[float(sp.re(ynew_clean)),float(sp.im(ynew_clean))],
 'new_observed_response':[float(sp.re(ynew)),float(sp.im(ynew))],
 'true_material_residual_upper':float(rn/2),
 'wrong_class_trace_lower':float(exact[1]['trace']),
 'wrong_class_residual_squared_margin':float(wrong_margin),
 'wrong_class_residual_squared_margin_rational':str(wrong_margin),
 'wrong_class_distance_lower_from_trace':float(sp.sqrt(exact[1]['trace'])),
 'baseline_inside_first_box':baseline_inside,'exact_ambiguous_alternative_inside_second_box':alternative_inside,
 'scoring_contract':'Only existing measured Yold, box endpoints, G,S,Bdict and fixed weight bank used before choosing the angle. Truth new response computed only after selection. Best angle means best among this finite menu of candidate weights/directions, not global optimal illumination.',
 'proof_contract':'Positive K by exact leading principal minors; old threshold negative and conditional receiver coefficient positive by exact rational signs; disk norm radii enclosed via integer-square-root rational brackets; clean prediction uses source-homogeneity q<=0; gap minus 2r_new positive; wrong class excluded by exact tr(H_aug)>r_joint^2 with nu=0; truth retained by known exact residual r_new/2.',
}
(folder/'active_ambiguity_closed_loop.json').write_text(json.dumps(result,indent=2))
np.savez(folder/'active_ambiguity_closed_loop_operators.npz',G=G,S=S,Bd=Bd)
print(json.dumps({k:v for k,v in result.items() if not k.endswith('_rational')},indent=2))
