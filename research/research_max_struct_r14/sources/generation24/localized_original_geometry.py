def ellipse(cx,cy,rx,ry,angle=0.):
    co,si=np.cos(angle),np.sin(angle);x=co*(xx-cx)+si*(yy-cy);y=-si*(xx-cx)+co*(yy-cy)
    return (x/rx)**2+(y/ry)**2<=1


def segment(ax,ay,bx,by,width):
    dx,dy=bx-ax,by-ay;t=np.clip(((xx-ax)*dx+(yy-ay)*dy)/(dx*dx+dy*dy),0,1)
    return (xx-ax-t*dx)**2+(yy-ay-t*dy)**2<=width*width


def material(case):
    rng=np.random.default_rng(case['material_seed']);f=case['family_index'];v=case['family_variant']
    H=float(rng.uniform(12,18));L=float(rng.uniform(2,5))
    uH=(H-12.)/6.;uL=(L-2.)/3.;regime=int(case['contrast_regime'])
    high,low=((4.+uH,2.+uL) if regime==0 else (15.+3.*uH,10.+2.*uL) if regime==1 else (H,L))
    angle=float(rng.uniform(-np.pi,np.pi));shift=rng.uniform(-.05,.05,2)
    h=high;m=np.full(xx.shape,.1);detail={'high':high,'low':low,'rotation_radians':angle,'shift':shift.tolist()}
    def put(mask,val):m[mask]=val
    if f==0:
        count=int(rng.integers(3,7));objs=[]
        for j in range(count):
            cx,cy=rng.uniform(-.29,.29,2);rx,ry=rng.uniform(.055,.13,2);a=float(rng.uniform(-np.pi,np.pi));val=h if j%2==0 else low
            put(ellipse(cx,cy,rx,ry,a),val);objs.append([cx,cy,rx,ry,a,val])
        detail['ellipses']=objs
    elif f==1:
        z=np.zeros(xx.shape);bumps=[]
        for j in range(8):
            cx,cy=rng.uniform(-.3,.3,2);width=float(rng.uniform(.07,.19));amp=float(rng.uniform(-1,1));z+=amp*np.exp(-((xx-cx)**2+(yy-cy)**2)/(2*width*width));bumps.append([cx,cy,width,amp])
        z=(z-z.min())/(z.max()-z.min());taper=np.clip((.40-np.sqrt((xx-shift[0])**2+(yy-shift[1])**2))/.09,0,1)
        m=.1+taper*(low+(h-low)*z-.1);detail['bumps']=bumps
    elif f==2:
        centers=rng.uniform(-.32,.32,(7,2));d=(xx[:,:,None]-centers[:,0])**2+(yy[:,:,None]-centers[:,1])**2;label=np.argmin(d,axis=2)
        vals=rng.uniform(low,h,7);vals[::3]=h;mask=ellipse(*shift,.37,.29,angle)
        put(mask,vals[label][mask]);detail.update(centers=centers.tolist(),levels=vals.tolist())
    elif f==3:
        co,si=np.cos(angle),np.sin(angle);x=co*(xx-shift[0])+si*(yy-shift[1]);y=-si*(xx-shift[0])+co*(yy-shift[1]);body=(abs(x)<.34)&(abs(y)<.27)
        put(body,h);holes=[]
        for cx,cy,r in [(-.15,-.1,.065),(.12,.11,.085),(.16,-.12,.045)]:
            cx+=float(rng.uniform(-.025,.025));cy+=float(rng.uniform(-.025,.025));hole=(x-cx)**2+(y-cy)**2<r*r;put(hole,.1);holes.append([cx,cy,r])
        put(body&(x<-.25),low);detail['holes_in_body_coordinates']=holes
    elif f==4:
        co,si=np.cos(angle),np.sin(angle);x=co*(xx-shift[0])+si*(yy-shift[1]);y=-si*(xx-shift[0])+co*(yy-shift[1]);body=(abs(x)<.34)&(abs(y)<.29)
        # A C-shaped rectangular notch and a separate closed oval hole.
        notch=(x>-.10)&(abs(y)<rng.uniform(.065,.10));closed=((x+.21)/.055)**2+((y-.14)/.06)**2<1
        put(body&~notch&~closed,h);put(body&(y<-.19),low)
    elif f==5:
        pts=np.array([[-.27,-.21],[.25,-.18],[-.03,.26]])+shift
        widths=[float(rng.uniform(.018,.033)),float(rng.uniform(.025,.045))]
        for j in [0,1]:put(segment(*pts[j],*pts[2],widths[j]),low if j==0 else h)
        for j,(cx,cy) in enumerate(pts):put(ellipse(cx,cy,float(rng.uniform(.07,.105)),float(rng.uniform(.055,.095)),angle+j),h if j!=1 else low)
        detail.update(junctions=pts.tolist(),bridge_halfwidths=widths)
    elif f==6:
        put(ellipse(*shift,.35,.285,angle),h)
        center=shift+rng.uniform(-.06,.06,2);put(ellipse(*center,.25,.18,angle+.3),low)
        inner=center+rng.uniform(-.035,.035,2);put(ellipse(*inner,.13,.09,angle-.4),h*.65)
        put(ellipse(*(inner+np.array([.02,-.01])),.055,.04,angle),.1)
        detail.update(layer_centers=[shift.tolist(),center.tolist(),inner.tolist()])
    else:
        centers=np.array([[-.24,-.20],[.19,.20],[-.21,.21],[.22,-.21]])+rng.uniform(-.035,.035,(4,2))
        for j,(cx,cy) in enumerate(centers):
            rx,ry=(.12,.085) if j<2 else(.05,.075);mask=ellipse(cx,cy,rx,ry,angle+j*.4)
            val=low if j%2 else h
            # One smooth internal gradient gives a mixed piecewise/continuous object.
            if j==1:put(mask,(val+2*(xx-cx))[mask])
            else:put(mask,val)
        detail['island_centers']=centers.tolist()
    return np.clip(m,.02,20),detail


def average(a,side):
    q=216//side;return a.reshape(side,q,side,q).mean(axis=(1,3)).ravel()


def acquisition(side,kind):
    p=grids[side];si=np.arange(0,12,2) if kind=='full6' else np.arange(12);ri=np.arange(24) if kind=='half12' else np.arange(48)
    return dict(G=p['G'],E0=p['E0'][:,si],S=p['S'][ri],coords=p['centers'],centers=p['centers'],sources=p['sources'][si],receivers=p['receivers'][ri],grid_shape=np.array([side,side]),k=float(p['k']),radius=float(p['radius']),cell_radius=float(p['radius'])),si,ri

