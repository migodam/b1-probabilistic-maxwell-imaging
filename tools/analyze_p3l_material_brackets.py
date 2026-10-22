"""Push saved projector-bracket diagnostics into the actual material metric.

Zero Maxwell actions. A nonzero coordinate bracket alone is not a material
obstruction when the decoder has a numerical kernel. No truth is loaded.
"""
from pathlib import Path
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from p3l.core import OUT, Book, write, save


def main():
    book = Book('C', 'material_bracket_saved_arrays', max_wall=300)
    status, error = 'COMPLETE', None
    rows = []
    try:
        folder = OUT / 'C/material_brackets'
        if folder.exists():
            raise FileExistsError('saved-material-bracket analysis already exists')
        folder.mkdir(parents=True)
        write(folder/'PROTOCOL.json', dict(
            parents=['DEV_000', 'DEV_028'], steps=[.01,.005,.0025],
            inputs='C/jobs/geometry_projector_DEV000_DEV028; all arrays already paid',
            material_metric='stack(Re Dg_z, Im Dg_z)/24; complex RMS on 576 cells',
            frame='eight-dimensional saved posterior-material-risk span',
            kernel_thresholds=[1e-8,1e-10,1e-12],
            proposition='For a C2 immersion g, pushforward Lie bracket is Dg times coordinate bracket: symmetric D2g terms cancel. At a numerical rank defect only the represented material diagnostic is asserted.',
            no_exact_fiber_claim=True, no_global_integrability_claim=True,
            GT_read=False, new_Maxwell_RHS=0))
        for parent in ['DEV_000','DEV_028']:
            base=OUT/'C/jobs/geometry_projector_DEV000_DEV028'/parent
            with np.load(base/'BASE_RAW.npz',allow_pickle=False) as z:
                Dc=z['material_derivative'][:,:64]
            D=np.vstack((Dc.real,Dc.imag))/24
            with np.load(base/'BASE_GEOMETRY.npz',allow_pickle=False) as z:
                V=z['Vw'].copy()
            sv=np.linalg.svd(D,compute_uv=False)
            U,ws,_=np.linalg.svd(D@V,full_matrices=False)
            weak_rank=int(np.sum(ws>ws[0]*1e-10))
            Q=U[:,:weak_rank]
            parent_rows=[]
            normals={}
            for h in [.01,.005,.0025]:
                with np.load(base/f'h{h:g}_DIRECTIONAL_PROJECTOR.npz',allow_pickle=False) as z:
                    rec={}
                    raw={}
                    for kind in ['actual','frozen']:
                        b=z[kind+'_bracket'].copy()
                        pushed=D@b
                        normal=pushed-Q@(Q.T@pushed)
                        tangent=Q@(Q.T@pushed)
                        rec[kind]=dict(material_bracket_RMS=float(np.linalg.norm(pushed)),
                            material_normal_RMS=float(np.linalg.norm(normal)),
                            material_tangent_RMS=float(np.linalg.norm(tangent)),
                            normal_fraction=float(np.linalg.norm(normal)/max(1e-300,np.linalg.norm(pushed))),
                            coordinate_bracket_norm=float(np.linalg.norm(b)),
                            material_speed_seed1=float(np.linalg.norm(D@V[:,0])),
                            material_speed_seed2=float(np.linalg.norm(D@V[:,1])))
                        raw[kind+'_material_bracket']=pushed
                        raw[kind+'_material_normal']=normal
                        normals[(kind,h)]=normal
                entry=dict(parent=parent,h=h,**rec)
                parent_rows.append(entry); rows.append(entry)
                save(folder/f'{parent}_h{h:g}.npz',**raw)
            differences=[]
            for a,b in [(.01,.005),(.005,.0025)]:
                differences.append(dict(coarse_h=a,fine_h=b,
                    actual_material_normal_difference_RMS=float(np.linalg.norm(normals['actual',a]-normals['actual',b])),
                    frozen_material_normal_difference_RMS=float(np.linalg.norm(normals['frozen',a]-normals['frozen',b]))))
            doc=dict(parent=parent,material_Dg_singular_values=sv,
                rank_by_relative_threshold={str(t):int(np.sum(sv>sv[0]*t)) for t in [1e-8,1e-10,1e-12]},
                material_weak_span_rank=weak_rank,material_weak_span_singular_values=ws,
                rows=parent_rows,successive_step_differences=differences,
                status='ACTUAL_MATERIAL_METRIC_FINITE_DIFFERENCE_DIAGNOSTIC',
                interpretation='Nonzero normal leakage concerns this chosen smooth projector distribution and tested pair. It does not prove absence of other charts or of finite-noise posterior regions.')
            write(folder/f'{parent}.json',doc)
        write(folder/'SUMMARY.json',dict(status='COMPLETE_ZERO_NEW_PHYSICS',rows=rows,
            original_geometry_job='geometry_projector_DEV000_DEV028',
            Maxwell_RHS=0,states=0,no_GT=True,
            restrictions='Two anchors, one pair each, finite-difference and floating-point evidence; no theorem of global nonintegrability or exact fiber.'))
        print(json.dumps([dict(parent=r['parent'],h=r['h'],normal_RMS=r['actual']['material_normal_RMS']) for r in rows]))
    except BaseException as exc:
        status,error='ERROR',repr(exc)
        raise
    finally:
        book.close(status,error)


if __name__=='__main__':
    main()
