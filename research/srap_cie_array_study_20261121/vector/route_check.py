from projection import *
mod=Model(12,'plane_grid');Y=t(np.load(ROOT/'data/plane_grid_00.npz')['Y']);val,U=torch.linalg.eigh(mod.S@mod.S.mH);sel=val>.01*val[-1];V=mod.S.mH@U[:,sel]/val[sel].sqrt()[None,:]
p=Projector(mod,V,Y,float(torch.linalg.norm(Y))*.1)
print(p.info,flush=True)
