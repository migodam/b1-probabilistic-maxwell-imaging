"""Frozen G_GRID_0125_V1 geometry; no observation or true-mask dependency."""
from itertools import product
import numpy as np

GRID_ID = 'G_GRID_0125_V1'
ROI = (-1.5, 1.5, 0.5, 1.25)
CELL = 0.125
SHAPE = (6, 24)  # row=y, column=x


def expand_support(support):
    a = np.asarray(support)
    if a.shape not in ((2, 2), (3, 3), SHAPE):
        raise ValueError('support must have shape (2,2), (3,3), or (6,24)')
    if not np.all((a == 0) | (a == 1)):
        raise ValueError('support must be binary')
    return np.repeat(np.repeat(a.astype(bool), SHAPE[0] // a.shape[0], 0), SHAPE[1] // a.shape[1], 1)


def _components(mask):
    todo = set(map(tuple, np.argwhere(mask)))
    groups = []
    while todo:
        seed = todo.pop()
        group, stack = {seed}, [seed]
        while stack:
            y, x = stack.pop()
            for q in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)):
                if q in todo:
                    todo.remove(q)
                    group.add(q)
                    stack.append(q)
        groups.append(group)
    return groups


def _material(value, support, name):
    a = np.asarray(value, dtype=float)
    if a.ndim == 0:
        return np.full(SHAPE, float(a))
    if a.shape == SHAPE:
        return a.copy()
    if a.shape == np.asarray(support).shape:
        return np.repeat(np.repeat(a, SHAPE[0] // a.shape[0], 0), SHAPE[1] // a.shape[1], 1)
    raise ValueError(name + ' must be scalar or match support/grid')


def disk_wall_gap(mask, center, radius):
    """Exact Euclidean distance from a closed disk to occupied closed rectangles."""
    m = expand_support(mask)
    p = np.asarray(center, dtype=float)
    if p.shape != (2,) or not np.all(np.isfinite(p)) or not np.isfinite(radius) or radius < 0:
        raise ValueError('invalid disk')
    y, x = np.where(m)
    if len(x) == 0:
        return float('inf')
    lo_x, lo_y = ROI[0] + x*CELL, ROI[2] + y*CELL
    dx = np.maximum(np.maximum(lo_x-p[0], p[0]-(lo_x+CELL)), 0)
    dy = np.maximum(np.maximum(lo_y-p[1], p[1]-(lo_y+CELL)), 0)
    return float(max(0, np.min(np.hypot(dx, dy))-radius))


def validate_geometry(support, epsinf=2.0, ell=0.0, target_center=None, target_radius=None):
    m = expand_support(support)
    e, l = _material(epsinf,support,'epsinf'), _material(ell,support,'ell')
    reasons = []
    occupied = _components(m)
    if len(occupied) > 1:
        reasons.append('occupied_support_disconnected')
    padded = np.pad(m, 1, constant_values=False)
    a,b,c,d = padded[:-1,:-1],padded[:-1,1:],padded[1:,:-1],padded[1:,1:]
    if np.any((a & d & ~b & ~c) | (b & c & ~a & ~d)):
        reasons.append('checkerboard_diagonal_contact')
    air = _components(~padded)
    holes = sum(not any(y in (0,padded.shape[0]-1) or x in (0,padded.shape[1]-1) for y,x in g) for g in air)
    area = float(m.sum()*CELL**2)
    if area > 2.25 + 1e-12:
        reasons.append('material_area_exceeded')
    if not np.all(np.isfinite(e[m])) or np.any((e[m]<1.5)|(e[m]>6)):
        reasons.append('epsinf_out_of_range')
    if not np.all(np.isfinite(l[m])) or np.any((l[m]<0)|(l[m]>0.25)):
        reasons.append('ell_out_of_range')
    # Air parameters are fixed by the model, regardless of scalar material input.
    e[~m], l[~m] = 1.0, 0.0
    gap = None
    if target_center is not None or target_radius is not None:
        p = np.asarray(target_center,dtype=float)
        r = target_radius
        if p.shape != (2,) or not np.all(np.isfinite(p)) or r is None or not np.isfinite(r):
            reasons.append('invalid_target')
        else:
            if not (-1<=p[0]<=1 and 1.6<=p[1]<=2.6 and 0.08<=r<=0.18):
                reasons.append('target_out_of_range')
            if r >= 0:
                gap = disk_wall_gap(m,p,r)
                if gap < 0.15-1e-12:
                    reasons.append('wall_target_gap_below_0.15')
    return dict(legal=not reasons,reasons=reasons,components=len(occupied),air_components=len(air),holes=int(holes),area=area,gap=gap,support=m,epsinf=e,ell=l,grid_id=GRID_ID)


def enumerate_supports(size=2, probability=0.65):
    """List (coarse_support, conditional_probability), including empty support."""
    if size not in (2,3) or not (0<probability<1):
        raise ValueError('size must be 2 or 3 and probability strictly between 0 and 1')
    legal = []
    for bits in product((False,True), repeat=size*size):
        support = np.asarray(bits,dtype=bool).reshape(size,size)
        if validate_geometry(support)['legal']:
            n = int(support.sum())
            legal.append((support,probability**n*(1-probability)**(size*size-n)))
    z = sum(p for _,p in legal)
    return [(s,p/z) for s,p in legal]


validate = validate_geometry
