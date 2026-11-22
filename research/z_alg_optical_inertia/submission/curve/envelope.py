"""Fourier lower envelope for arbitrary two-by-two real and optical weights."""

from itertools import product
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh
from scipy.fft import fft2
from flint import acb, acb_mat, arb

from native import ball, gershgorin, high, low, mid, star
from search import h, times


class FourierEnvelope:
    def __init__(self, geometry):
        self.g = geometry
        self.roots = [acb(arb(2)*k/geometry.pad).exp_pi_i() for k in range(geometry.pad)]
        self.soph = {}
        self.dirichlet = {}
        self.ss = geometry.sensor_ball*star(geometry.sensor_ball)
        self.precise_sensor_blocks = None
        nblk = geometry.block_side
        indices=np.arange(geometry.side//nblk)
        f = np.exp(2j*np.pi*np.arange(geometry.pad)[:,None]*indices[None,:]/geometry.pad).sum(axis=1)
        self.numerical_dirichlet = np.stack([
            f*np.exp(2j*np.pi*np.arange(geometry.pad)*len(indices)*j/geometry.pad) for j in range(nblk)])
        n, p = geometry.side, geometry.pad
        self.numerical_ports = np.stack([
            fft2((geometry.sensor_star*(geometry.blocks==a)[:,None]).reshape(n,n,-1),
                 s=(p,p),axes=(0,1))/p for a in range(geometry.block_count)])
        self.numerical_ss_blocks = np.stack([
            geometry.sensor[:,geometry.blocks==a] @ geometry.sensor_star[geometry.blocks==a,:]
            for a in range(geometry.block_count)])
        self.numerical_grams = {}

    def precise_fourier_ports(self, modes):
        """Return F* S* restricted to each block, caching unique Fourier entries."""
        g, roots = self.g,self.roots
        needed = [z for z in modes if z not in self.soph]
        if needed:
            a = sorted({j for j,k in needed})
            b = sorted({k for j,k in needed})
            accum = {z:[[] for j in range(g.block_count)] for z in needed}
            for block in range(g.block_count):
                end = g.side//g.block_side
                row, col = divmod(block,g.block_side)
                rows = range(row*end,(row+1)*end)
                cols = range(col*end,(col+1)*end)
                yf = acb_mat(len(a),end,[roots[-k*i % g.pad] for k in a for i in rows])
                xf = acb_mat(end,len(b),[roots[-k*i % g.pad] for i in cols for k in b])
                for rec in range(len(g.sensor)):
                    view = acb_mat(end,end,[g.sensor_ball[rec,i*g.side+j].conjugate()
                                            for i in rows for j in cols])
                    transformed=yf*view*xf/g.pad
                    for j,k in needed:
                        accum[j,k][block].append(transformed[a.index(j),b.index(k)])
            self.soph.update(accum)
        return [self.soph[z] for z in modes]

    def precise_entry(self, dy, dx):
        """Return the four blockwise F*F entries for a Fourier difference."""
        key = (dy,dx)
        if key not in self.dirichlet:
            nblk=self.g.block_side
            n,p=self.g.side//nblk,self.g.pad
            y = [sum((self.roots[dy*k%p] for k in range(j*n,(j+1)*n)),acb(0)) for j in range(nblk)]
            x = [sum((self.roots[dx*k%p] for k in range(j*n,(j+1)*n)),acb(0)) for j in range(nblk)]
            self.dirichlet[key] = [yy*xx/(p*p) for yy in y for xx in x]
        return self.dirichlet[key]

    def weighted_sensor_gram(self, alpha, gap):
        matrix = acb_mat(len(self.g.sensor),len(self.g.sensor))
        if self.precise_sensor_blocks is None:
            self.precise_sensor_blocks = []
            nrec = len(self.g.sensor)
            for j in range(self.g.block_count):
                pixels = np.flatnonzero(self.g.blocks==j).tolist()
                sub = acb_mat(nrec,len(pixels),[self.g.sensor_ball[i,k]
                                               for i in range(nrec) for k in pixels])
                self.precise_sensor_blocks.append(sub*star(sub))
        for a,d in zip(self.precise_sensor_blocks,alpha,strict=True):
            matrix += a/(d-gap)
        return matrix

    @staticmethod
    def numeric_coeffs(w, a, beta, shift):
        normalized = np.asarray(w.real)/w.scale
        t = np.mean(w.optical)/w.scale
        z = ((normalized-1)*w.cap + 1j*(np.asarray(w.optical)/w.scale-t))/2
        diag = normalized-1-2*shift*z.real
        alpha = a + diag-beta*abs(z)**2
        return t,z,alpha

    def numerical_matrix(self,w,cut,beta,shift,gap):
        g=self.g
        t,z,diagonal=self.numeric_coeffs(w,0.,beta,shift)
        iso = w.is_global and np.allclose(np.array(w.real),w.scale,rtol=0.,atol=0.)
        symbol = 1-w.cap*g.spec.real+t*g.spec.imag
        cross = g.spec-shift
        score = symbol if iso else symbol-abs(cross)**2/beta
        if cut not in self.numerical_grams:
            mode1 = np.arange(-cut,cut+1) % g.pad
            y,x = np.meshgrid(mode1,mode1,indexing='ij')
            iy,ix = y.ravel(),x.ravel()
            dy=(iy[None,:]-iy[:,None])%g.pad
            dx=(ix[None,:]-ix[:,None])%g.pad
            blocks=np.array([self.numerical_dirichlet[j//g.block_side,dy]
                             *self.numerical_dirichlet[j%g.block_side,dx]/g.pad**2
                             for j in range(g.block_count)])
            ports=self.numerical_ports[:,iy,ix,:]
            omitted=np.ones(g.spec.shape,dtype=bool)
            omitted[iy,ix]=False
            self.numerical_grams[cut] = iy,ix,blocks,ports,omitted
        iy,ix,blocks,ports,omitted = self.numerical_grams[cut]
        a = float(np.floor(128*score[omitted].min())/128)
        alpha=a+diagonal
        if np.min(alpha) <= gap+1.e-10:
            return None
        v, score1, f = symbol[iy,ix]-a,score[iy,ix]-a,cross[iy,ix]
        if min(abs(score1)) < 1.e-11:
            return None
        r = len(iy)
        rec = len(g.sensor)
        invd=1/(alpha-gap)
        src = np.einsum('a,aij->ij',invd,ports)
        gg = np.einsum('a,aij->ij',invd,blocks)
        ss = np.einsum('a,aij->ij',invd,self.numerical_ss_blocks)
        scale = w.scale
        if iso:
            invmat = np.block([[gg,src],[src.conj().T,ss]])
            idx = np.arange(r)
            invmat[idx,idx] += 1/score1
        else:
            gz = np.einsum('a,aij->ij',invd*z,blocks)
            gzz = np.einsum('a,aij->ij',invd*abs(z)**2,blocks)
            pz = np.einsum('a,aij->ij',invd*z.conj(),ports)
            invmat = np.block([[gg,gz,src],[gz.conj().T,gzz,pz],
                               [src.conj().T,pz.conj().T,ss]])
            idx=np.arange(r)
            det=beta*score1
            invmat[idx,idx] += beta/det
            invmat[idx+r,idx+r] += v/det
            invmat[idx,idx+r] += f.conj()/det
            invmat[idx+r,idx] += f/det
        invmat[-rec:,-rec:] += scale*np.eye(rec)
        return invmat,int(np.count_nonzero(score1<0)),tuple(zip(iy.tolist(),ix.tolist())),iso,a

    def select(self,w,nu):
        """Select a numerically well-separated proof proposal from nested Fourier squares."""
        for cut in (4,6,8,10,12):
            for gap in (.03125,.015625,.0078125,.001953125):
                values = ((1.,0.),) if w.is_global else (
                    (4.,.25),(2.,.25),(8.,.25),(16.,.25),(1.,.25),(4.,.375),(8.,.375))
                for beta,shift in values:
                    result = self.numerical_matrix(w,cut,beta,shift,gap)
                    if result is None:
                        continue
                    matrix,n,ret,iso,a=result
                    need = max(0,n-nu)
                    if need == 0:
                        return a,beta,shift,gap,ret,iso,need
                    val=eigh(h(matrix),eigvals_only=True,check_finite=False,subset_by_index=(need-1,need-1))
                    if val[0] < -1.e-7:
                        return a,beta,shift,gap,ret,iso,need
        raise ArithmeticError('Fourier lower envelope did not locate the native positive gap.')

    def verify(self,w,nu):
        start,cpu=perf_counter(),process_time()
        envelope, beta0, shift0, gap0, retained, iso, need = self.select(w,nu)
        g=self.g
        a,beta,shift,gap,scale=map(arb,(envelope,beta0,shift0,gap0,w.scale))
        cap=arb(w.cap)
        t=sum((arb(l) for l in w.optical),arb(0))/(g.block_count*scale)
        wr=[arb(x)/scale for x in w.real]
        z=[acb((r-1)*cap,arb(l)/scale-t)/2 for r,l in zip(wr,w.optical,strict=True)]
        if iso:
            wr=[arb(1)]*g.block_count
            z=[acb(0)]*g.block_count
        alpha=[a+r-1-2*shift*d.real-beta*(d.real*d.real+d.imag*d.imag)
               for r,d in zip(wr,z,strict=True)]
        min_alpha=min(value.lower() for value in alpha)
        if not min_alpha>gap or (iso and not all(d==0 for d in z)):
            raise ArithmeticError('Invalid identity floor for native Fourier proof.')
        keep=set(retained)
        scores={}
        minimum=arb(1000000000)
        for j,srow in enumerate(g.spec_ball):
            for k,g0 in enumerate(srow):
                w0=1-cap*g0.real+t*g0.imag-a
                f=g0-shift
                score=w0 if iso else w0-(f.real*f.real+f.imag*f.imag)/beta
                if (j,k) in keep:
                    if not (score<0 or score>0):
                        raise ArithmeticError('Retained Fourier block determinant is unresolved.')
                    scores[j,k]=(w0, f, score)
                    continue
                minimum=min(minimum,score.lower())
        if not minimum>=0:
            raise ArithmeticError('An omitted Fourier block is not positive semidefinite.')
        n_neg = sum(score < 0 for _,_,score in scores.values())
        if need != max(0,n_neg-nu):
            raise ArithmeticError('The compressed inertia dimension disagrees with exact retained signs.')
        rec=len(g.sensor)
        r=len(retained)
        off=r if iso else 2*r
        dim=off+rec
        if need == 0:
            return dict(positive_gap_lower=gap0,upper_count=nu,proof_dimension=dim,
                        retained_fourier_modes=r,omitted_symbol_slack_lower=low(minimum),
                        lower_identity_floor=low(min_alpha),**times(start,cpu))
        matrix=acb_mat(dim,dim)
        receiver=self.precise_fourier_ports(retained)
        gram={}
        z2=[d.real*d.real+d.imag*d.imag for d in z]
        vinv=[1/(a0-gap) for a0 in alpha]
        for i,(iy,ix) in enumerate(retained):
            for j,(jy,jx) in enumerate(retained):
                delta=(jy-iy)%g.pad,(jx-ix)%g.pad
                if delta not in gram:
                    rates=[x*d for x,d in zip(self.precise_entry(*delta),vinv,strict=True)]
                    gram[delta]=(
                        sum(rates,acb(0)),
                        sum((x*d for x,d in zip(rates,z,strict=True)),acb(0)),
                        sum((x*d.conjugate() for x,d in zip(rates,z,strict=True)),acb(0)),
                        sum((x*d for x,d in zip(rates,z2,strict=True)),acb(0)))
                x00,x01,x10,x11=gram[delta]
                matrix[i,j]=x00
                if not iso:
                    matrix[i,j+r],matrix[i+r,j],matrix[i+r,j+r]=x01,x10,x11
            for j in range(rec):
                blocks=[x[j]/(d-gap) for x,d in zip(receiver[i],alpha,strict=True)]
                row=sum(blocks,acb(0))
                matrix[i,off+j],matrix[off+j,i]=row,row.conjugate()
                if not iso:
                    row=sum((x*d.conjugate() for x,d in zip(blocks,z,strict=True)),acb(0))
                    matrix[i+r,off+j],matrix[off+j,i+r]=row,row.conjugate()
            v, f, score=scores[iy,ix]
            if iso:
                matrix[i,i]+=1/score
            else:
                det=beta*score
                matrix[i,i]+=beta/det
                matrix[i+r,i+r]+=v/det
                matrix[i,i+r]+=f.conjugate()/det
                matrix[i+r,i]+=f/det
        sensor_inverse_gram=self.ss/(alpha[0]-gap) if iso else self.weighted_sensor_gram(alpha,gap)
        for i in range(rec):
            for j in range(rec):
                matrix[off+i,off+j]=sensor_inverse_gram[i,j]
            matrix[off+i,off+i]+=scale
        _,basis=eigh(h(mid(matrix)),check_finite=False,subset_by_index=(0,need-1))
        basis=ball(basis)
        upper=gershgorin(star(basis)*(matrix*basis),upper=True)
        if not upper<0:
            raise ArithmeticError('Outward compressed inertia count unresolved.')
        return dict(positive_gap_lower=gap0,upper_count=nu,proof_dimension=dim,
                    retained_fourier_modes=r,omitted_symbol_slack_lower=low(minimum),
                    lower_identity_floor=low(min_alpha),proof_negative_subspace_dimension=need,
                    proof_negative_subspace_upper=high(upper),
                    symbol_floor=envelope,young_beta=beta0,symbol_center=shift0,
                    **times(start,cpu))
