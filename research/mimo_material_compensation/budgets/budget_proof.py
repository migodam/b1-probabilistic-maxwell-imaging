"""A posteriori bounds for the fixed nominal released multi-port operator."""
from quick import np, lu_factor, lu_solve
from verified_wave import VerifiedWave, fro, opnorm, rd, up


def down(value):
    return float(np.nextafter(float(value) - rd(3)*abs(value), -np.inf))


def norm_lower(a):
    a=np.asarray(a)
    return max(0., down(float(np.linalg.norm(a.ravel()))/(1.+rd(2*a.size))))


class VerifiedState:
    def __init__(self, wave, rec, exc, x, full_k=False):
        self.x, self.rec, self.exc = x, rec, exc
        active=np.flatnonzero(x)
        d=len(active)
        self.active=active
        n=len(x)
        self.wave=wave
        if d==0:
            self.j=np.zeros_like(exc)
            self.p=np.zeros((len(rec),exc.shape[1]),complex)
            self.u,self.e=rec,exc
            self.u_error=self.e_error=self.response_error=0.
            self.k_total=wave.gnorm
            self.log=dict(active_dimension=0, k_total_upper=self.k_total,
                          response_error_upper=0., inverse_defect_upper=0.)
            return
        dense = np.arange(n)
        gre,er=self.green_block(active,dense)
        gle,el=self.green_block(dense,active)
        g=gre[:,active]
        _,eg=self.green_block(active,active)
        cap=float(x.max())
        amount=x[active]
        matrix=np.eye(d)-amount[:,None]*g
        arith=rd(3)*(np.sqrt(d)+cap*fro(g))
        matrix_error=up(cap*eg+arith)
        factor=lu_factor(matrix,check_finite=False)
        inverse=lu_solve(factor,np.eye(d),check_finite=False)
        inverse_norm=opnorm(inverse)
        defect=up(fro(np.eye(d)-matrix@inverse)+matrix_error*inverse_norm+
                  rd(d+4)*(np.sqrt(d)+fro(matrix)*fro(inverse)))
        if defect>=1:
            raise ArithmeticError("Failed active-state inverse verification")
        inv_bound=up(inverse_norm/(1.-defect))
        inv_error=up(inv_bound*defect)
        weighted=inverse*amount
        w_error=up(cap*inv_error+rd(2)*fro(weighted))
        rhs=amount[:,None]*exc[active]
        j=lu_solve(factor,rhs,check_finite=False)
        j_error=up(inv_bound*(fro(matrix@j-rhs)+matrix_error*fro(j)
                    +rd(d+4)*(fro(matrix)*fro(j)+cap*fro(exc[active]))))
        self.j=np.zeros_like(exc)
        self.j[active]=j
        self.p=rec[:,active]@j
        self.response_error=up(opnorm(rec[:,active])*j_error+rd(d)*fro(rec[:,active])*fro(j))
        self.e=exc+gle@j
        self.e_error=up((opnorm(gle)+el)*j_error+el*fro(j)+rd(d+4)*(fro(exc)+fro(gle)*fro(j)))
        c_by_w=rec[:,active]@weighted
        cw_error=up(w_error*fro(rec[:,active])+rd(d)*fro(rec[:,active])*fro(weighted))
        self.u=rec+c_by_w@gre
        self.u_error=up(cw_error*(opnorm(gre)+er)+opnorm(c_by_w)*er+
                        rd(d+4)*(fro(rec)+fro(c_by_w)*fro(gre)))
        self.cross_left,self.cross_right,self.weighted=gle,gre,weighted
        self.cross_left_error,self.cross_right_error,self.weighted_error=el,er,w_error
        self.k_total=self.k_upper(np.ones(n,bool), np.ones(n,bool), dense=full_k)
        self.log=dict(active_dimension=d, inverse_defect_upper=defect, inverse_norm_upper=inv_bound,
                      response_error_upper=self.response_error, field_e_error_upper=self.e_error,
                      receive_u_error_upper=self.u_error, k_total_upper=self.k_total)

    def green_block(self, rows, cols):
        wave=self.wave
        yr,xr=divmod(rows,wave.side)
        yc,xc=divmod(cols,wave.side)
        out=wave.lag[(yr[:,None]-yc[None,:])%wave.p,(xr[:,None]-xc[None,:])%wave.p]
        left,right=wave.L[rows],wave.R[:,cols]
        lr=left@right
        error=up(np.sqrt(len(rows)*len(cols))*wave.lagerr+
                 rd(left.shape[1])*fro(left)*fro(right)+rd(2)*(fro(out)+fro(lr)))
        return out+lr,error

    def k_upper(self, rows, cols, dense):
        if not len(self.active):
            return self.wave.gnorm
        left,right=self.cross_left[rows],self.cross_right[:,cols]
        el,er,ew=self.cross_left_error,self.cross_right_error,self.weighted_error
        if not dense:
            return up(self.wave.gnorm+(opnorm(left)+el)*(opnorm(self.weighted)+ew)*(opnorm(right)+er))
        direct,ed=self.green_block(np.flatnonzero(rows),np.flatnonzero(cols))
        d=len(self.active)
        wright=self.weighted@right
        ewr=up(ew*(opnorm(right)+er)+opnorm(self.weighted)*er+rd(d)*fro(self.weighted)*fro(right))
        correction=left@wright
        ec=up(el*(opnorm(wright)+ewr)+opnorm(left)*ewr+rd(d)*fro(left)*fro(wright))
        return up(opnorm(direct+correction)+ed+ec+rd(2)*(fro(direct)+fro(correction)))

    def jacobian(self, layout):
        tx,rx=layout["tx"],layout["rx"]
        u,e=self.u[rx],self.e[:,tx]
        product=(u[:,None,:]*e.T[None,:,:]).reshape(len(rx)*len(tx),-1)
        j=np.concatenate((product.real,product.imag),axis=0)
        error=up(self.u_error*(fro(e)+self.e_error)+fro(u)*self.e_error+rd(4)*fro(u)*fro(e))
        return j, error, u, e


class VerifiedNativeState(VerifiedState):
    """Matrix-free audit of a dense material on the native grid."""
    def __init__(self,wave,rec,exc,x,proposal):
        self.x,self.rec,self.exc,self.wave=x,rec,exc,wave
        self.active=np.flatnonzero(x)
        inv,invlog=wave.inverse(x)
        self.j=proposal["j"]
        err_j=wave.solve_error(x,self.j,x[:,None]*exc,up(rd(3)*x.max()*fro(exc)),inv)
        self.p=rec@self.j
        self.response_error=up(rd(len(x))*fro(rec)*fro(self.j)+opnorm(rec)*err_j)
        g=wave.apply(self.j,strict=True)
        self.e=exc+g
        self.e_error=up(wave.full_green_err*fro(self.j)+wave.gnorm*err_j+rd(3)*(fro(exc)+fro(g)))
        if "u" in proposal:
            self.u=proposal["u"]
            # For A* U*=C*, compute the transpose of the strict transform
            # without symmetry assumptions on the stored port coefficients.
            v=x[:,None]*self.u.conj().T
            a=np.zeros((v.shape[1],wave.p,wave.p),complex)
            a[:,:wave.side,:wave.side]=v.T.reshape((-1,wave.side,wave.side))
            tr=wave.transform(a)*wave.sym.conj()
            gp=wave.transform(tr,True)[:,:wave.side,:wave.side].reshape(v.shape[1],-1).T
            gp+=wave.R.conj().T@(wave.L.conj().T@v)
            residual=self.u.conj().T-gp-rec.conj().T
            err=up(wave.full_green_err*fro(v)+wave.gnorm*rd(3)*fro(v)+
                   rd(4)*(fro(self.u)+fro(gp)+fro(rec))+fro(residual))
            self.u_error=up(inv*err)
        self.k_total=up(wave.gnorm*inv)
        self.log=dict(inverse=invlog,active_dimension=len(self.active),
                      response_error_upper=self.response_error,
                      field_e_error_upper=self.e_error,
                      receive_u_error_upper=getattr(self,"u_error",None), k_total_upper=self.k_total)

    def k_upper(self,rows,cols,dense):
        return self.k_total
