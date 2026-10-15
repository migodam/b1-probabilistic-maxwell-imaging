"""Exact radial volume integration; no truth available to online callers."""
import numpy as np
from b1_joint.material import decode,material_distance
from .target import to_x

def metrics(a,b):
    da,db=decode(to_x(a)),decode(to_x(b));edges=sorted(set(da['radii'])|set(db['radii'])|{0.})
    er=ei=ec=inter=union=void_inter=void_union=roi_ec=0.
    for lo,hi in zip(edges[:-1],edges[1:]):
        mid=(lo+hi)/2;ea=da['eps_layers'][np.searchsorted(da['radii'],mid)];eb=db['eps_layers'][np.searchsorted(db['radii'],mid)];v=hi**3-lo**3
        er+=v*abs(ea.real-eb.real);ei+=v*abs(ea.imag-eb.imag);ec+=v*abs(ea-eb)
        va=abs(ea-1)>1e-12;vb=abs(eb-1)>1e-12
        inter+=v*(va and vb);union+=v*(va or vb);void_inter+=v*((not va) and (not vb));void_union+=v*((not va) or (not vb))
        roi_ec+=v*abs(ea-eb)*((not va) or (not vb))
    # Declared fixed material dynamic ranges, vacuum included.
    norm=.5*(er/1.8+ei/.16)
    return dict(real_mae=er,imag_mae=ei,complex_mae=ec,normalized_mae=norm,complex_nmae=ec/abs(1.8+.16j),support_iou=inter/union if union else 1.,cavity_iou=void_inter/void_union if void_union else 1.,central_radius_error=abs(a[0]-b[0]),ring_halfwidth_error=abs(a[1]-b[1]),boundary_distance=max(abs(a[0]-b[0]),abs(a[1]-b[1])),material_l2=material_distance(to_x(a),to_x(b)),cavity_region_complex_mae=roi_ec/void_union if void_union else None)

def distance(a,b):return material_distance(to_x(a),to_x(b))

def structural_difference(a,b):return abs(a[0]-b[0])>=.05 or abs(a[1]-b[1])>=.025

def nested_select(candidates,maximum=8):
    if not candidates:return []
    # Legal shared base is always first; otherwise lowest residual replaces it.
    selected=[0];remaining=list(range(1,len(candidates)))
    while remaining and len(selected)<maximum:
        i=max(remaining,key=lambda j:min(distance(candidates[j]['p'],candidates[k]['p']) for k in selected))
        selected.append(i);remaining.remove(i)
    return selected
