"""Generator-only geometry; never import this module in an inverse worker."""
import numpy as np

FAMILIES = ('annulus','slotted_shell','multi_island','curved_bars','comb',
            'coated_inclusion','perforated_star','bridge')


def draw_material(rng, family, regime, high_chi, side=24):
    """One independent draw, no forward-based rejection or reroll.

    Coordinates are fractions of the ROI halfwidth. Every feature width is
    fixed before generation; the 24-grid rendering itself is the finite truth.
    """
    if family not in FAMILIES or regime not in ('low','high','mixed'):
        raise ValueError('Unknown frozen generator family/regime')
    axis = (np.arange(side)+.5)*2/side-1
    yy,xx = np.meshgrid(axis,axis,indexing='ij')
    tx,ty = rng.uniform(-.08,.08,2)
    angle = rng.uniform(0,2*np.pi)
    sx,sy = rng.uniform(.9,1.1,2)
    u = ((xx-tx)*np.cos(angle)+(yy-ty)*np.sin(angle))/sx
    v = (-(xx-tx)*np.sin(angle)+(yy-ty)*np.cos(angle))/sy
    r=np.hypot(u,v);theta=np.arctan2(v,u)
    parameters=dict(translation=[float(tx),float(ty)],angle=float(angle),
                    scales=[float(sx),float(sy)])
    if family == 'annulus':
        outer=rng.uniform(.56,.69);inner=rng.uniform(.27,.36)
        support=(r<outer)&(r>inner)
        parameters.update(outer=float(outer),inner=float(inner))
    elif family == 'slotted_shell':
        outer=rng.uniform(.60,.70);inner=rng.uniform(.34,.40)
        slot=rng.uniform(.12,.20)
        support=(r<outer)&(r>inner)&~((u>0)&(abs(v)<slot))
        parameters.update(outer=float(outer),inner=float(inner),slot_halfwidth=float(slot))
    elif family == 'multi_island':
        sep=rng.uniform(.30,.39);rad=rng.uniform(.19,.25)
        support=((u+sep)**2+(v+.08)**2<rad**2)|((u-sep)**2+(v-.08)**2<(rad*.86)**2)
        support|=((u+.02)**2+(v-.48)**2<(rad*.65)**2)
        parameters.update(separation=float(sep),radius=float(rad))
    elif family == 'curved_bars':
        rad=rng.uniform(.46,.56);width=rng.uniform(.13,.17)
        support=((abs(np.hypot(u+.14,v)-rad)<width)&(u>-.10)&(v<.51))
        support|=((abs(np.hypot(u-.17,v+.08)-rad*.72)<width*.76)&(u<.12)&(v>-.47))
        parameters.update(radius=float(rad),halfwidth=float(width))
    elif family == 'comb':
        shift=rng.uniform(-.035,.035);width=rng.uniform(.10,.13)
        support=(abs(v+.48)<width)&(abs(u)<.65)
        for center in [-.50,-.17,.17,.50]:
            support|=(abs(u-center-shift)<width*.72)&(v>-.47)&(v<.52)
        parameters.update(shift=float(shift),halfwidth=float(width))
    elif family == 'coated_inclusion':
        a=rng.uniform(.56,.68);b=rng.uniform(.40,.49)
        support=(u/a)**2+(v/b)**2<1
        parameters.update(major=float(a),minor=float(b))
    elif family == 'perforated_star':
        phase=rng.uniform(-.2,.2);outer=rng.uniform(.51,.60)
        support=r<outer*(1+.22*np.cos(5*(theta+phase)))
        support&=~(((u-.16)**2+(v+.05)**2<.16**2)|((u+.22)**2+(v-.08)**2<.12**2))
        parameters.update(phase=float(phase),radius=float(outer))
    else:  # bridge: two thick bodies connected by a narrow neck
        rad=rng.uniform(.26,.32);sep=rng.uniform(.37,.43);neck=rng.uniform(.10,.14)
        support=((u-sep)**2+v**2<rad**2)|((u+sep)**2+v**2<rad**2)
        support|=(abs(u)<sep)&(abs(v)<neck)
        parameters.update(radius=float(rad),separation=float(sep),neck_halfwidth=float(neck))
    high_region = support & (u+.23*v>0)
    if family == 'coated_inclusion':
        high_region=support&(((u-.09)/.27)**2+((v+.04)/.21)**2<1)
    m=np.full((side,side),.1)
    if regime=='low':m[support]=2.5
    elif regime=='high':m[support]=high_chi
    else:
        m[support]=2.5;m[high_region]=high_chi
    if not support.any() or (regime=='mixed' and
            (not high_region.any() or not (support&~high_region).any())):
        raise ValueError('Frozen draw is degenerate; fail without reroll')
    return m.ravel(),parameters
