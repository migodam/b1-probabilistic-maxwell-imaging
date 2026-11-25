"""Shared definitions for the prospective limited-aperture scalar-TM experiment."""
import json
import sys
from hashlib import blake2b
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor, lu_solve, eigh
from scipy.special import hankel1, jv

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY / "ip_guided_v4/boxes"))
from run import Geometry, coordinates, native_geometry

# The registered limited-aperture menu uses the first half of the 32-way dictionary.
RX96 = np.arange(16)
TX96 = np.arange(16)
OLD = 0
NOISE_DRAW = 202611260464
GEOMETRY_SEED = 202611260443
BACKGROUND_WIDTH = 0.0015
FOREGROUND_WIDTH = 0.001
N_TEST = 20


def atomic_write(path, payload):
    part = path.with_suffix('.part')
    part.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    part.replace(path)


def save(path, **arrays):
    part = path.with_suffix('.part.npz')
    np.savez_compressed(part, **arrays)
    part.replace(path)


def require_absent(*paths):
    for path in paths:
        if path.exists():
            raise FileExistsError(path)


def digest(array):
    array = np.ascontiguousarray(array)
    return blake2b(array.tobytes(), digest_size=20).hexdigest()


def clocks():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])


def herm(matrix):
    return (matrix+matrix.conj().T)/2


def geometry(side):
    full = native_geometry()
    if side == 96:
        return Geometry(full.S[RX96], full.B[:,TX96], side, spectrum=full.spectrum)
    if side != 32:
        raise ValueError(side)
    with np.load(STUDY/'synthetic/public/coarse_green.npz') as f:
        native = f['G']
    xy = coordinates(side)
    angles = np.arange(32)*np.pi/16
    direction = np.column_stack((np.cos(angles), np.sin(angles)))
    wave = 2*np.pi*400e6/3e8
    b = np.exp(1j*wave*(xy@direction[TX96].T))*np.sqrt(10/32)
    radius = 2/side/np.sqrt(np.pi)
    pos = 3.75*direction[RX96]
    d = np.linalg.norm(pos[:,None] - xy[None,:], axis=2)
    s = np.sqrt(20/32)*.5j*np.pi*wave*radius*jv(1,wave*radius)*hankel1(0,wave*d)
    return Geometry(s,b,side,G=native)


def forward(g, s, b, chi, *, auxiliary=False):
    a = np.eye(len(chi))-chi[:,None]*g
    lu = lu_factor(a, check_finite=False)
    current = lu_solve(lu, chi[:,None]*b, check_finite=False)
    measured = s@current
    if not auxiliary:
        return measured
    field = b+g@current
    receive = lu_solve(lu, s.T, trans=1, check_finite=False).T
    born = s@(chi[:,None]*b)
    return measured, field, receive, current, born, lu


def energy_gram(geom):
    denominator = np.vdot(geom.B[:,0],geom.B[:,0]).real
    return herm(geom.B.conj().T@geom.B/denominator)


def normalize(vec, gram, total=1.):
    return vec*np.sqrt(total/np.vdot(vec,gram@vec).real)


def shapes(side, index, kind):
    """Independent deterministic shape/pose parameters, assigned before any data."""
    rng = np.random.default_rng(np.random.SeedSequence([GEOMETRY_SEED,index]))
    cx,cy = rng.uniform(-.18,.18,size=2)
    angle = rng.uniform(-np.pi,np.pi)
    peak = rng.uniform(3.2,5.5)
    r = rng.uniform(.315,.4)
    rx,ry = rng.uniform(.89,1.12,size=2)*r
    xy = coordinates(side)
    x0,y0=xy.T
    x = (x0-cx)*np.cos(angle)+(y0-cy)*np.sin(angle)
    y = -(x0-cx)*np.sin(angle)+(y0-cy)*np.cos(angle)
    d = (x/rx)**2+(y/ry)**2
    if kind=='solid':
        fg = d <= 1
    elif kind=='ring':
        fg = (d<=1)&(d>=.4**2)
    elif kind=='double':
        fg = ((x+.62*r)/(.55*rx))**2+(y/(.73*ry))**2<=1
        fg |= ((x-.58*r)/(.53*rx))**2+((y-.12*r)/(.78*ry))**2<=1
    elif kind=='austria':
        fg = (d<=.88**2)&(d>=.49**2)&~((x>0)&(y>.22*r)&(y<.72*r))
        fg |= ((x-.96*r)/(.24*rx))**2+((y-.54*r)/(.26*ry))**2<=1
        fg |= ((x-.94*r)/(.21*rx))**2+((y+.53*r)/(.24*ry))**2<=1
    bg = (~fg)&(((x+1.18*r)**2+(y-.9*r)**2<(.25*r)**2)|
                  ((x-.94*r)**2+(y+.93*r)**2<(.32*r)**2))
    active = np.flatnonzero(fg|bg)
    q = np.maximum(0.,1-d)
    chi = np.where(fg,peak*(.85+.15*q),.045)
    return active, fg[active], bg[active], chi[active], np.column_stack([x[active]/r,y[active]/r])


def profiles(xy, foreground):
    x,y = xy.T
    centers = np.array([[xx,yy] for xx in np.linspace(-.75,.75,4) for yy in np.linspace(-.75,.75,4)])
    dist2 = ((xy[:,None]-centers[None,:])**2).sum(axis=-1)
    funcs = np.column_stack((np.ones(len(x)),x,y,np.exp(-dist2/.095)))
    funcs *= foreground[:,None]
    return funcs[:,np.linalg.norm(funcs,axis=0)>1.e-8]


def pair_jacobian(field, receive, funcs):
    return np.einsum('pi,is,ij->psj',receive,field,funcs,optimize=True)
