#!/usr/bin/env python3
"""Offline certificate of true-material feasibility in the exact saved Fourier model."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import argparse
import gc
import json
from pathlib import Path
import resource
import sys

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.linalg.blas import dgemm
from flint import acb, acb_mat, arb, ctx

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
DEFAULT_OUT = STUDY/'results/theory/fine_feasibility.json'
sys.path.insert(0, str(STUDY))
from core.certificate import ball, frobenius_squared, safe_float
from core.inertia import clock, timing

ARRAYS = ('PW180_10x20', 'PW360_10x20')
GUARD = arb(1)/10**6
BITS = 128
CHUNK = 128
EPS = arb(2)**-53
UNDERFLOW_GUARD = arb(2)**-1022


def hi(v):
    return safe_float(v,upper=True)


def lo(v):
    return safe_float(v,upper=False)


def gamma(k):
    return arb(k)*EPS/(1-arb(k)*EPS)


def binary_frob_ceil(x):
    """Upper norm via exact component maxima of a binary64 array."""
    re = arb(float(np.abs(x.real).max()))
    im = arb(float(np.abs(x.imag).max()))
    return (arb(x.size)*(re*re+im*im)).sqrt().upper()


def acb_value(z):
    return acb(float(z.real),float(z.imag))


class ExactFourierGeometry:
    def __init__(self, path):
        with np.load(path,allow_pickle=False) as data:
            spectral = data['spectrum']
            apertures = [(data[f'{name}__S'],data[f'{name}__B']) for name in ARRAYS]
        self.pad = spectral.shape[0]
        assert spectral.shape==(self.pad,self.pad)
        self.side = int(np.sqrt(apertures[0][1].shape[0]))
        assert self.side**2==apertures[0][1].shape[0] and self.pad>=2*self.side-1
        # NumPy-compatible inverse DFT along x and y, each normalized by pad.
        horizontal = [acb.dft([acb_value(z) for z in row],inverse=True) for row in spectral]
        vertical = [acb.dft([row[j] for row in horizontal],inverse=True) for j in range(self.pad)]
        self.exact = [vertical[j][i] for i in range(self.pad) for j in range(self.pad)]
        self.centers = np.fromiter((complex(z) for z in self.exact),complex,len(self.exact))
        self.error_max = max(abs(z-acb_value(w)).upper() for z,w in zip(self.exact,self.centers,strict=True))
        self.kernel_center_max = max(abs(acb_value(w)).upper() for w in self.centers)
        self.S, ports = np.unique(np.concatenate([a[0] for a in apertures],axis=0),axis=0,return_inverse=True)
        b, sources = np.unique(np.concatenate([a[1].T for a in apertures],axis=0),axis=0,return_inverse=True)
        self.B = b.T
        num_rx = np.cumsum([a[0].shape[0] for a in apertures])[:-1]
        num_tx = np.cumsum([a[1].shape[1] for a in apertures])[:-1]
        self.inds = dict(zip(ARRAYS,zip(np.split(ports,num_rx),np.split(sources,num_tx),strict=True),strict=True))

    def active_block_indices(self, active):
        y, x = np.divmod(active,self.side)
        return ((y[:,None]-y[None,:]) % self.pad)*self.pad+((x[:,None]-x[None,:]) % self.pad)


class StableActiveModel:
    def __init__(self, geometry, material):
        self.geo = geometry
        self.active = np.flatnonzero(material)
        self.m = material[self.active]
        if not np.min(self.m)>0 or np.max(self.m)>4:
            raise ValueError('Active cells must have 0 < χ <= 4.')
        self.n = len(self.active)
        self.I = geometry.active_block_indices(self.active)
        self.S = geometry.S[:,self.active]
        self.B = geometry.B[self.active]
        self.A = -self.m[:,None]*geometry.centers[self.I]
        self.A.flat[::self.n+1] += 1.
        self.lu = lu_factor(self.A,check_finite=False)

    def inverse_bound(self):
        started=clock()
        n=self.n
        approx=lu_solve(self.lu,np.eye(n,dtype=complex),check_finite=False)
        Rr,Ri=approx.real,approx.imag
        Ar,Ai=self.A.real,self.A.imag
        real_residual=dgemm(-1.,Rr,Ar)+dgemm(1.,Ri,Ai)
        real_residual.flat[::n+1]+=1.
        imag_residual=dgemm(-1.,Rr,Ai)+dgemm(-1.,Ri,Ar)
        approximate_defect=binary_frob_ceil(real_residual+1j*imag_residual)
        # Four independent real GEMMs plus final subtraction/addition. The
        # envelope uses the standard binary64 summation/product error model
        # with more operations than a full dot-product dependency chain.
        R_norm,A_norm=map(binary_frob_ceil,(approx,self.A))
        product_guard=gamma(4*n+32)*(4*R_norm*A_norm+arb(n).sqrt())
        product_guard+=arb((4*n+32)*n*n)*UNDERFLOW_GUARD
        # Coefficient enclosure for the exact A versus binary64 assembly A0.
        # Bound real and imaginary perturbations jointly, including forming I-MG.
        m_max=arb(float(self.m.max()))
        formation_guard = gamma(8)*(arb(2).sqrt()+2*m_max*self.geo.kernel_center_max)+16*UNDERFLOW_GUARD
        enclosure_frob=n*(formation_guard+m_max*self.geo.error_max)
        total=(approximate_defect+product_guard+R_norm*enclosure_frob).upper()
        if not total<1:
            raise ArithmeticError('Inverse Neumann test failed on the exact Fourier operator.')
        upper=(R_norm/(1-total)).upper()
        return upper,dict(approximate_inverse_norm_upper=hi(R_norm),
            binary64_inverse_defect_upper=hi(approximate_defect), product_rounding_guard_upper=hi(product_guard),
            A_exact_minus_binary64_frobenius_upper=hi(enclosure_frob),
            inverse_neumann_defect_upper=hi(total), inverse_spectral_norm_upper=hi(upper),
            state_singular_value_min_lower=lo(1/upper), **timing(started))

    def enclose_field(self, inv_bound, chunk):
        started=clock()
        n=self.n
        state=lu_solve(self.lu,self.m[:,None]*self.B,check_finite=False)
        adj=lu_solve(self.lu,self.S.T,trans=1,check_finite=False).T
        X,W,SB,BB=map(ball,(state,adj,self.S,self.B))
        ms=[arb(float(v)) for v in self.m]
        state_defect_entries=[]
        receiver_defect=SB-W
        for lo_i in range(0,n,chunk):
            end_i=min(n,lo_i+chunk)
            rows=range(lo_i,end_i)
            block=self.I[lo_i:end_i]
            G=acb_mat(len(rows),n,[self.geo.exact[int(k)] for k in block.ravel()])
            GX=G*X
            state_defect_entries.extend(ms[i]*(BB[i,s]+GX[i-lo_i,s])-X[i,s]
                                        for i in rows for s in range(BB.ncols()))
            small_adj=acb_mat(W.nrows(),len(rows),[W[j,i]*ms[i] for j in range(W.nrows()) for i in rows])
            receiver_defect += small_adj*G
        defect=acb_mat(n,BB.ncols(),state_defect_entries)
        primal=frobenius_squared(defect).sqrt().upper()
        dual=frobenius_squared(receiver_defect).sqrt().upper()
        error=(dual*inv_bound*primal).upper()
        corrected=SB*X+W*defect
        return corrected,error,dict(active_state_columns=BB.ncols(), used_receiver_rows=W.nrows(),
            primal_defect_frobenius_upper=hi(primal), receiver_adjoint_defect_frobenius_upper=hi(dual),
            uncorrected_state_error_frobenius_upper=hi(inv_bound*primal),
            corrected_port_remainder_frobenius_upper=hi(error),**timing(started))


def norm(v):
    return frobenius_squared(v).sqrt()


def audit_observation(geometry, corrected, global_error, observation):
    file=observation['file']
    with np.load(HERE/'public'/file,allow_pickle=False) as raw:
        Y=ball(raw['Y'])
    with np.load(HERE/'evaluation'/file,allow_pickle=False) as raw:
        saved=ball(raw['Y_fine_clean'])
    rows,cols=geometry.inds[observation['geometry']]
    port=acb_mat(len(rows),len(cols),[corrected[int(i),int(j)] for i in rows for j in cols])
    exact_y_distance=norm(Y-port)
    left=(exact_y_distance.lower()-global_error).lower()
    right=(exact_y_distance.upper()+global_error).upper()
    solver_delta=(norm(saved-port).upper()+global_error).upper()
    y_norm=norm(Y)
    radius=arb(float(observation['radius_noise_fro']))
    # Round the enlarged declared value upward to binary64 for future proof replay.
    padded_float=hi(radius+GUARD*y_norm)
    padded=arb(padded_float)
    needed=max(arb(0),((right-radius)/y_norm.lower()).upper())
    return dict(id=observation['id'],geometry=observation['geometry'],noise_fraction=observation['noise_fraction'],
        norm_Y_upper=hi(y_norm),sample_radius=hi(radius),sample_radius_binary64=float(observation['radius_noise_fro']),
        true_forward_to_Y_lower=max(0.,lo(left)), true_forward_to_Y_upper=hi(right),
        true_forward_to_saved_numerical_clean_upper=hi(solver_delta),
        true_forward_to_saved_numerical_clean_relative_upper=hi(solver_delta/y_norm.lower()),
        raw_status='inside' if right<=radius else 'outside' if left>radius else 'unresolved',
        absorb_solver_relative_allowance=1.e-6, minimum_sufficient_relative_guard_upper=hi(needed),
        inflated_radius=padded_float, inside_inflated_radius=bool(right<=padded),
        inflated_radius_margin_lower=lo(padded-right))


def proof_margin_audit(observations):
    file=STUDY/'results/theory/fine_verified.json'
    proofs=json.loads(file.read_text(encoding='utf8'))['cases']
    by_id={d['id']:d for d in observations}
    rows=[]
    for entry in proofs:
        if entry['id'] not in by_id:
            continue
        info=by_id[entry['id']]
        if float(entry['noise_radius']) != info['sample_radius_binary64']:
            raise ArithmeticError('Native fine certificate used a different noise ball.')
        radius=info['inflated_radius']
        margin=arb(float(entry['distance_lower']))-arb(radius)
        rows.append(dict(id=entry['id'], tested_cap=entry['cap'], profile=entry['profile'],
            distance_lower=entry['distance_lower'], inflated_radius=radius,
            true_material_admitted=info['inside_inflated_radius'], noise_only_margin_lower=entry['margin_lower'],
            inflated_margin_lower=lo(margin),
            excludes_on_inflated_radius=bool(margin>0 and info['inside_inflated_radius'])))
    return rows


def write(path, val):
    path.parent.mkdir(parents=True,exist_ok=True)
    pending=path.with_suffix('.partial')
    pending.write_text(json.dumps(val,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    pending.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--chunk',type=int,default=CHUNK)
    parser.add_argument('--object',action='append')
    args=parser.parse_args()
    ctx.prec,ctx.threads=BITS,1
    start=clock()
    t=clock()
    geometry=ExactFourierGeometry(HERE/'public/fine_geometry.npz')
    notes=json.loads((HERE/'public/manifest.json').read_text(encoding='utf8'))['observations']
    groups={k:[v for v in notes if v['object_id']==k] for k in sorted({v['object_id'] for v in notes})}
    report=dict(geometry='synthetic/public/fine_geometry.npz', binary64_exact_data=True,
        material_labelling='offline after native fine-grid data-only certificates; exact saved discrete operator',
        inverse_roundoff_model='IEEE binary64 four real DGEMM dot-product bound, Arb enclosed coefficient perturbation',
        fourier_operator='top-left 96² embedding and restriction with stored exact binary64 192² Fourier multiplier',
        precision_bits=BITS, arb_threads=1, row_chunk=args.chunk,
        max_ifft_coefficient_error_from_binary64_center=hi(geometry.error_max),
        shared_setup=timing(t), objects={}, observations=[])
    for object_id, obs in groups.items():
        if args.object and object_id not in args.object:
            continue
        t=clock()
        with np.load(HERE/'evaluation'/f'{object_id}.npz',allow_pickle=False) as src:
            m=src['m_fine']
        model=StableActiveModel(geometry,m)
        est,inverse=model.inverse_bound()
        corrected,bound,forward=model.enclose_field(est,args.chunk)
        data=[audit_observation(geometry,corrected,bound,o) for o in obs]
        report['objects'][object_id]={'nonzero_true_fine_cells':model.n,
            'inverse_bound':inverse,'adjoint_correction':forward,**timing(t)}
        report['observations'].extend(data)
        report['cost']={**timing(start),'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss /
                        (1024**2 if sys.platform=='darwin' else 1024)}
        report['checked_certificate_margins']=proof_margin_audit(report['observations'])
        write(args.output,report)
        print(object_id,model.n,'raw',[d['raw_status'] for d in data],
            'all inflated',all(d['inside_inflated_radius'] for d in data),
            'max solver delta',max(d['true_forward_to_saved_numerical_clean_upper'] for d in data),flush=True)
        del model, corrected
        gc.collect()


if __name__=='__main__':
    main()
