"""Public, frozen R2 wall and target families. No observation or scoring input."""
import numpy as np
from .physics import WallGrid,Target
def axes(eta=.125):
    step=eta/2
    return np.linspace(-1.5,1.5,int(round(3/step))+1),np.linspace(.5,1.25,int(round(.75/step))+1)
def micro_wall(parent,eta,texture):
    x,y=axes(eta);xx,yy=np.meshgrid((x[:-1]+x[1:])/2,(y[:-1]+y[1:])/2)
    occupied=(xx>=-.75)&(xx<.75)
    if parent==1:
        occupied &= ~((yy>=1)&(xx<-.5))
        occupied &= ~((yy>=.75)&(yy<1)&(xx>=-.25)&(xx<0))
    cx=-1.5+(np.floor((xx+1.5)/.25)+.5)*.25
    cy=.5+(np.floor((yy-.5)/.25)+.5)*.25
    mean=3.35+.25*np.cos(np.pi*cx/.75)+.15*(cy-.875)/.375
    loss=.065+.01*np.cos(2*np.pi*cx)
    ix=np.floor((xx+1.5)/(eta/2)+1e-7).astype(int)
    iy=np.floor((yy-.5)/(eta/2)+1e-7).astype(int)
    if texture==0:p=(-1.)**(ix+iy)
    elif texture==1:p=(-1.)**(ix+iy//2)
    elif texture==-1:p=np.zeros_like(xx)
    else:raise ValueError('texture')
    eps=np.where(occupied,mean+.8*p,1.);ell=np.where(occupied,loss,0.)
    return WallGrid(x,y,eps,ell)
def public_wall(g):
    x=np.linspace(-1.5,1.5,13);y=np.linspace(.5,1.25,4)
    xx,yy=np.meshgrid((x[:-1]+x[1:])/2,(y[:-1]+y[1:])/2)
    support=(abs(xx)<.75)
    if g//2==1:
        support=(abs(xx)<.75)&~((yy>=1)&(xx<-.5))&~((yy>=.75)&(yy<1)&(xx>=-.25)&(xx<0))
    if g//2==2:support=np.ones_like(xx,dtype=bool)
    return WallGrid(x,y,np.where(support,[3.0,3.8][g%2],1.),np.where(support,.065,0.))
def target_defs():
    train=[dict(pose=[-.4,1.8],radius=.08,contrast=1.,loss=0.,shape='disk',orientation=0.),
           dict(pose=[.5,2.4],radius=.16,contrast=6.,loss=.1,shape='ellipse',orientation=float(np.pi/2))]
    held=[]
    for s in range(2):
        for j,(px,py) in enumerate([(-.65,1.75),(.15,1.75),(.65,1.75),(-.65,2.45),(.15,2.45),(.65,2.45)]):
            held.append(dict(pose=[px,py],radius=[.08,.12,.16][(j+s)%3],contrast=[1.,3.,6.][(2*j+s)%3],loss=[0.,.1,.2][(j+2*s)%3],shape=['disk','ellipse'][s],orientation=float((j%2)*np.pi/2)))
    return train,held
def make_target(d):return Target(**d)
def white(z):
    a=np.asarray(z).ravel();return np.r_[a.real,a.imag]/.01
def whiteJ(z):return np.vstack((np.asarray(z).real,np.asarray(z).imag))/.01
