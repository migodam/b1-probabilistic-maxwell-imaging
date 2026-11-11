"""Explicit [0,8] adapter over existing R8 true-vector physical/PR primitives.

The shared R8 default [0,4] optimizer is not modified. This adapter supplies
bounds to the caller's optimization loop and adds counts for actual inherited
matrix actions. Physics, gradients, material objective, current PR and prior
are inherited verbatim.
"""
from research.research_max_struct_r8.src.vector_pr3d import VectorPR3D


class HighContrastVectorPR3D(VectorPR3D):
    material_bounds=(0.,8.)

    def __init__(self,ns=6,width=.6):
        super().__init__(ns=ns,width=width,lossy=False)
        # Existing R8 optional exact value cache: gradient calls retain their
        # paid fields, so an accepted-state audit need not repeat its LU.
        self.cache_physical_values=True

    def physical(self,z,Y=None,gradient=True):
        lu=self.counts['physical_LU'];adj=self.counts['adjoint_RHS']
        result=super().physical(z,Y,gradient)
        if self.counts['physical_LU']>lu:
            self.counts['physical_G_columns']=self.counts.get('physical_G_columns',0)+self.P
            self.counts['physical_S_columns']=self.counts.get('physical_S_columns',0)+self.P
        if self.counts['adjoint_RHS']>adj:
            self.counts['physical_S_adjoint_columns']=self.counts.get('physical_S_adjoint_columns',0)+self.P
        return result

    def material(self,z,J,kind,beta):
        result=super().material(z,J,kind,beta)
        self.counts['material_G_columns']=self.counts.get('material_G_columns',0)+self.P
        return result

    def prepare(self,Y):
        result=super().prepare(Y)
        self.counts['APCS_range_SVD']=self.counts.get('APCS_range_SVD',0)+(2 if self.whiten_apcs else 0)
        return result

    def current_pr(self,z,Y,metric,beta,mf,global_pr):
        result=super().current_pr(z,Y,metric,beta,mf,global_pr)
        cols=self.bases[mf].shape[1]
        for name,value in {'current_full_matrix_builds':1,'current_KB_columns':cols,
                           'current_SB_columns':cols,'current_KJd_columns':self.P,
                           'current_SJd_columns':self.P}.items():
            self.counts[name]=self.counts.get(name,0)+value
        return result
