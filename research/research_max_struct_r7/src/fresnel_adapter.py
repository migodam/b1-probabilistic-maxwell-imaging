"""2001 official Fresnel measured TM adapter; no pixel truth or generated data.

Unknown interior horn field is approximated by a cylindrical source normalized
from target-absent, on-axis opposite-receiver E only. All complex measurements
are conjugated from exp(+iwt) to the inherited outgoing H1 convention.
"""
from pathlib import Path
import numpy as np
from scipy.special import hankel1
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_z_alg_r4_diversity.src.physics import build_operators


def load_file(path,frequency=2.):
    a=np.loadtxt(path,comments='#');a=a[np.isclose(a[:,2],frequency)]
    total=a[:,3]+1j*a[:,4];inc=a[:,5]+1j*a[:,6]
    Y=np.full((72,36),np.nan+1j*np.nan);I=Y.copy()
    for row,y,e in zip(a,total-inc,inc):
        j,k=int(row[1])-1,int(row[0])-1
        if np.isfinite(Y[j,k]):raise ValueError('Duplicate measured key')
        Y[j,k]=np.conj(y);I[j,k]=np.conj(e)
    return a,Y,I


def prepare(path,frequency=2.,side=24,width=.16):
    a,Y,I=load_file(path,frequency)
    k=2*np.pi*frequency*1e9/299792458.
    h=width/side;axis=(np.arange(side)+.5)*h-width/2
    yy,xx=np.meshgrid(axis,axis,indexing='ij');centers=np.c_[xx.ravel(),yy.ravel()]
    sa=np.arange(36)*np.pi/18;ra=np.arange(72)*np.pi/36
    sources=.72*np.c_[np.cos(sa),np.sin(sa)];receivers=.76*np.c_[np.cos(ra),np.sin(ra)]
    public=build_operators(centers,sources,receivers,.45*h,k)
    opposite=(2*np.arange(36)+36)%72
    measured_axis=I[opposite,np.arange(36)]
    if not np.isfinite(measured_axis).all():raise ValueError('Missing incident calibration channel')
    # Source/receiver on-axis distance is de+dr; the horn field amplitude is
    # approximated using the same 2D outgoing cylindrical propagation.
    source_gain=measured_axis/(.25j*hankel1(0,k*(.72+.76)))
    public['E0']=public['E0']*source_gain[None,:]
    public['grid_shape']=np.array([side,side]);public['centers']=centers
    train_sources=np.arange(9)
    common=np.flatnonzero(np.isfinite(Y[:,train_sources]).all(axis=1))
    hold_receivers=common[np.arange(len(common))%5==0]
    train_receivers=np.setdiff1d(common,hold_receivers)
    train_mask=np.zeros(Y.shape,bool);train_mask[np.ix_(train_receivers,train_sources)]=True
    withheld=np.isfinite(Y)&~train_mask
    selected=dict(public,E0=public['E0'][:,train_sources],S=public['S'][train_receivers])
    model=MaxwellModel(selected)
    return dict(model=model,public=public,Y=Y,I=I,Y_train=Y[np.ix_(train_receivers,train_sources)],
                train_sources=train_sources,train_receivers=train_receivers,
                train_mask=train_mask,withheld_mask=withheld,centers=centers,
                source_gain=source_gain,common_receivers=common)
