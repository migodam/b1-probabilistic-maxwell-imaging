"""Fixed, truth-independent Green projections and absolute energy identities.

No network, fit, threshold selection, or physical solve is performed here.
Columns are independent source currents; real material vectors are embedded in
the same complex right-domain Green bases. SVD thresholds are frozen at .1 and
1e-12. Numerical intersections use the historical 1e-9 constraint tolerance.
"""
from pathlib import Path
import numpy as np
from scipy.linalg import svd, qr


def squared_norm(x):
    return float(np.vdot(x, x).real)


class FixedEnergyGeometry:
    """Reusable fixed projections; accepts (576,) or (576, sources)."""
    @classmethod
    def from_cache(cls, path, G=None, S=None, sensitivity=False):
        with np.load(Path(path), allow_pickle=False) as z:
            return cls(z['in_V'], z['in_s'], z['ob_V'], z['ob_s'], G, S, sensitivity)

    def __init__(self, in_V, in_s, ob_V, ob_s, G=None, S=None, sensitivity=False):
        self.G, self.S = G, S
        self.n = len(in_s)
        self.I = np.asarray(in_V)[:, np.asarray(in_s) >= in_s[0] * .1]
        self.Os = np.asarray(ob_V)[:, np.asarray(ob_s) >= ob_s[0] * .1]
        self.Ow = np.asarray(ob_V)[:, (np.asarray(ob_s) < ob_s[0] * .1) & (np.asarray(ob_s) > ob_s[0] * 1e-12)]
        self.O = np.concatenate((self.Os, self.Ow), axis=1)
        u, d, _ = svd(np.concatenate((self.I, self.O), axis=1), full_matrices=False, check_finite=False)
        self.B = u[:, d > 1e-9]
        # Numerical intersection I_weak ∩ O_weak, constraints in the 39-space.
        _, d, vh = svd(self.I.conj().T @ self.Ow, full_matrices=True, check_finite=False)
        self.IwOw = self.Ow @ vh[np.count_nonzero(d > 1e-9):].conj().T
        self.bands = {
            'P_in_strong19': (self.I, False), 'P_in_weak557': (self.I, True),
            'S_strong9': (self.Os, False), 'S_weak39': (self.Ow, False),
            'S_null528': (self.O, True), 'S_visible48': (self.O, False),
            'D509': (self.B, True), 'B67': (self.B, False),
        }
        self.metadata = dict(n=self.n, tau=.1, numerical_rank_rtol=1e-12,
            intersection_constraint_tolerance=1e-9,
            ranks={k:self.n-q.shape[1] if comp else q.shape[1] for k,(q,comp) in self.bands.items()},
            intersection_weak_weak_rank=self.IwOw.shape[1],
            D_definition='orthogonal complement of span(P_in strong19, S visible48); numerical weak/null intersection',
            material_scope='complex Green right-domain embedding of real material error, not the real material Jacobian')
        if self.metadata['ranks'] != dict(P_in_strong19=19,P_in_weak557=557,S_strong9=9,S_weak39=39,S_null528=528,S_visible48=48,D509=509,B67=67):
            raise ValueError('Frozen Green rank contract changed')
        self.sensitivity={}
        if sensitivity:
            for tau in (.03,.3):
                vi=np.asarray(in_V)[:,np.asarray(in_s)>=in_s[0]*tau]
                vo=np.asarray(ob_V)[:,np.asarray(ob_s)>=ob_s[0]*tau]
                u,d,_=svd(np.concatenate((vi,self.O),axis=1),full_matrices=False,check_finite=False)
                union=u[:,d>1e-9]
                self.sensitivity[f'tau{tau:g}']=(vi,vo,union)
            self.metadata['sensitivity']={k:dict(in_strong=v.shape[1],ob_strong=o.shape[1],dark=self.n-b.shape[1],bright=b.shape[1],training_changed=False)for k,(v,o,b)in self.sensitivity.items()}

    @staticmethod
    def _part(x, q, complement=False):
        y = q @ (q.conj().T @ x)
        return x-y if complement else y

    def project(self, x, band):
        return self._part(np.asarray(x), *self.bands[band])

    def energies(self, x, twofold=True):
        """Squared absolute norms; all twofold definitions stay separately named."""
        x = np.asarray(x)
        total = squared_norm(x)
        e = {'total':total}
        for name,(q,comp) in self.bands.items():
            v=squared_norm(q.conj().T @ x)
            e[name] = max(0.,total-v) if comp else v
        if self.G is not None: e['G_field'] = squared_norm(self.G @ x)
        if self.S is not None: e['S_field'] = squared_norm(self.S @ x)
        for tag,(vi,vo,union)in self.sensitivity.items():
            gi=squared_norm(vi.conj().T@x);so=squared_norm(vo.conj().T@x);bright=squared_norm(union.conj().T@x)
            e.update({tag+'__Gin_strong':gi,tag+'__Gin_weak':max(0.,total-gi),
                      tag+'__S_strong':so,tag+'__S_weak':max(0.,e['S_visible48']-so),
                      tag+'__dark':max(0.,total-bright),tag+'__bright':bright})
        if twofold:
            ip=[self._part(x,self.I), self._part(x,self.I,True)]
            op=[self._part(x,self.Os),self._part(x,self.Ow),self._part(x,self.O,True)]
            for a,ib in enumerate(('strong','weak')):
                for b,ob in enumerate(('strong','weak','null')):
                    e[f'I_to_O__{ib}__{ob}']=squared_norm(self._part(ip[a], [self.Os,self.Ow,self.O][b], b==2))
                    e[f'O_to_I__{ib}__{ob}']=squared_norm(self._part(op[b],self.I,a==1))
                    e[f'intersection__{ib}__{ob}']=0.
            e['intersection__weak__weak']=squared_norm(self.IwOw.conj().T @ x)
            e['intersection__weak__null']=e['D509']
            e['intersection_remainder']=max(0.,total-e['intersection__weak__weak']-e['intersection__weak__null'])
        return e

    def transition(self, before, after, twofold=True):
        """E_after-E_before = E_update + 2 Re <P before,P update>.

        The signed cross term is obtained by polarization, including sequential
        operators. This identity also applies to G/S output energies.
        """
        a=self.energies(before,twofold); b=self.energies(after,twofold)
        u=self.energies(np.asarray(after)-np.asarray(before),twofold)
        return {k:dict(before=a[k],after=b[k],delta=b[k]-a[k],update=u[k],
                       signed_cross=b[k]-a[k]-u[k]) for k in a}
