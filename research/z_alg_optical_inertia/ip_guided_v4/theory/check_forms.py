"""Verify the noncommuting background, optical sign, incidence congruence and singular control."""
from pathlib import Path
from fractions import Fraction
import json
from flint import acb, acb_mat, arb, ctx

def adj(a): return a.transpose().conjugate()
def eye(n): return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
def diag(values):return acb_mat([[v if i==j else 0 for j in range(len(values))] for i,v in enumerate(values)])
def block(a,b,c,d):return acb_mat([[z for x in row for z in x] for mats in [(a,b),(c,d)] for row in ([list(m[i,k] for k in range(m.ncols())) for m in mats] for i in range(mats[0].nrows()))])
def inv(a):return a.solve(eye(a.nrows()))
def im(a):return (a-adj(a))/acb(0,2)
def check(a):
    assert all(z.contains(0) for z in a.entries())
    bound=max(float(abs(z).upper()) for z in a.entries());assert bound<1e-30
    return bound
def main():
    ctx.prec=160
    g=acb_mat([[.5j,.25+.125j],[.25+.125j,.375j]])
    m=diag([2.,.5]);h=diag([.75,.25]);w=diag([.5,2.]);v=diag([-1.,.5]);a=h*w*h
    s=acb_mat([[.25+.5j,.5-.125j],[.25j,.625]])
    b=acb_mat([[1.,.5j],[.25-.5j,1.]]);y=acb_mat([[.125j,-1],[1-.25j,.5]])
    u=eye(2)-m*g; ec=inv(eye(2)-g*m)*b;fm=s*m*ec;yc=y-fm
    t=inv(eye(2)-g*m)*g;p=s*inv(u)
    k=adj(p)*p+w-adj(t)*a*t+im(v*t)
    l=adj(p)*yc+adj(t)*a*ec+acb(0,.5)*v*ec;c=adj(yc)*yc-adj(ec)*a*ec
    k0=adj(s)*s+adj(u)*w*u-adj(g)*a*g+im(v*g)
    l0=adj(s)*y+adj(u)*w*m*b+adj(g)*a*b+acb(0,.5)*v*b;c0=adj(y)*y+adj(b)*(m*w*m-a)*b
    tri=block(u,-m*b,acb_mat(2,2),eye(2))
    bordered=block(k,-l,-adj(l),c)
    checks={'bordered_congruence':check(adj(tri)*bordered*tri-block(k0,-l0,-adj(l0),c0))}
    x=diag([2.5,.375]);j=inv(eye(2)-x*g)*x*b;e=b+g*j;q=j-m*e
    checks['physical_center_order']=check(e-ec-t*q);checks['optical']=check(im(adj(q)*v*e))
    checks['receiver']=check(s*j-fm-p*q)
    loss=acb(1,.25);j=inv(eye(2)-loss*x*g)*(loss*x*b);e=b+g*j;q=j-loss*m*e
    z=adj(q)*v*e
    checks['rotated_loss_optical']=check(im(z)+(z+adj(z))/8)
    # A physical, singular K. Real regular contrast x=-1 attains the lower bound.
    fr=Fraction; H=fr(5,4);inc=fr(3,4);port=fr(3,4);obs=fr(-5,4)
    K=port**2+1-H**2; L=port*obs+H*inc
    C1=obs**2-inc**2; C2=fr(4);cur=-inc/(1+H)
    assert K==L==0 and C1==1 and (port*cur-obs)**2+C2==5
    out={'rigor':'160-bit complex balls for noncommuting two-cell identities; exact rational one-cell singular comparison','zero_identity_radii_upper':checks,'singular':{'state':0,'linear_term':0,'source_spectrum':[4,1],'integral_count_squared_distance':5,'squared_distance_at_x_minus_one':5}}
    Path(__file__).with_name('check_forms.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
