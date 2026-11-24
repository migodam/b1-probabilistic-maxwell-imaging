"""Data-only branch-box menu shared by the pruned and unpruned controls."""

import numpy as np

from native_box import Branch


def as_branch(name,caps):
    return Branch(name,tuple(tuple(float(x) for x in r) for r in np.atleast_2d(caps)))


def expand(branch,side):
    caps=np.asarray(branch.caps,float)
    if side%caps.shape[0] or side%caps.shape[1]:
        raise ValueError('Declared box tiles do not partition the solver grid.')
    return np.repeat(np.repeat(caps,side//caps.shape[0],axis=0),side//caps.shape[1],axis=1).ravel()


def tile_field(m,n):
    side=int(np.sqrt(len(m)))
    if side**2!=len(m) or side%n:
        raise ValueError('Prefix must be a square material grid divisible by the tile count.')
    k=side//n
    return m.reshape(n,k,n,k).mean(axis=(1,3))


def propose_branches(material,current=None,*,tile_side=4,tile_count=4,foreground=4.,background=.05,
                     wildcard=20.):
    """Prefix-to-box map. All hyperparameters are fixed before consuming the prefix."""
    m=np.maximum(np.asarray(material,float).ravel(),0.)
    response=tile_field(m*m,tile_side)
    if current is not None:
        energy=tile_field(np.mean(abs(np.asarray(current))**2,axis=1),tile_side)
        if energy.max()>0:
            response=response/max(response.max(),1e-30)+.15*energy/energy.max()
    rank=np.argsort(-response.ravel(),kind='stable')
    first=np.zeros(tile_side**2,bool)
    first[rank[:tile_count]]=True
    first=first.reshape(tile_side,tile_side)
    second=np.roll(first,tile_side//2,axis=1)
    if np.array_equal(first,second):
        second=np.roll(first,1,axis=0)
    if np.array_equal(first,second):
        second=np.roll(first,1,axis=1)
    if np.array_equal(first,second) or not 0<tile_count<tile_side**2:
        raise ValueError('Require two distinct nonempty localized hypotheses.')
    make=lambda name,mask,bg:as_branch(name,np.where(mask,foreground,bg))
    return {'branches': [as_branch('uniform_low_050',[[.5]]),
                         make('ranked_support',first,background),
                         make('shifted_support',second,background),
                         as_branch('wildcard',[[wildcard]])],
            'development_hard_controls':[make('ranked_hard_support',first,0.),make('shifted_hard_support',second,0.)],
            'tile_score':response.tolist(), 'tile_count':tile_count,'tile_side':tile_side,
            'foreground':foreground,'background':background}


def primary_menu(material,current=None):
    """The same four exact boxes are offered to both timing arms."""
    return propose_branches(material,current,tile_side=4,tile_count=1,foreground=4.,background=0.)


def central_safety_control():
    mask=np.zeros((8,8))
    mask[1:7,1:7]=4.
    return as_branch('offline_central_6_of_8_safety',mask)
