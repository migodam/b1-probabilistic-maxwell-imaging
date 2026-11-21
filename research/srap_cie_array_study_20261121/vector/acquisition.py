"""Report independently measurable port distances and polarization norms."""
from physics import ports,WIDTH,ROOT,write,np
out={}
for g in ['plane_grid','hemi']:
    out[g]={}
    for n,phase,role in [(8,.18,'transmit'),(16,.63,'receive')]:
        p,v=ports(n,g,phase);dist=np.linalg.norm(p,axis=1);gap=np.linalg.norm(np.maximum(abs(p)-WIDTH/2,0),axis=1)
        out[g][role]=dict(positions=p.tolist(),polarizations=v.tolist(),center_range=[dist.min(),dist.max()],cube_boundary_range=[gap.min(),gap.max()],outer_sphere_range=[dist.min()-np.sqrt(3)*WIDTH/2,dist.max()-np.sqrt(3)*WIDTH/2],max_orthonormality_error=float(abs(v@v.transpose(0,2,1)-np.eye(2)).max()))
write(ROOT/'acquisition.json',out)
