"""One formula-selected full-wave interaction graph, then exact paired rays.

No nearest-distance/gradient/GT threshold. Each cell picks exactly one partner
maximizing the nonnegative squared mixed second-order response. At most576
unordered pairs survive the union. Pause after this single dictionary.
"""
from pathlib import Path
import sys,time
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
from paired_material_rays import run
from finite_coordinate_paths import fixed_edges

def select(model,C,O,Q,BO,BQ,public,cost,directory):
    start=time.perf_counter();cpu=time.process_time()
    bo=np.diag(BO).real;bq=np.diag(BQ).real
    norm=abs(C)**2*bo[:,None]*bq[None,:]+abs(C.T)**2*bo[None,:]*bq[:,None]
    norm+=2*np.real(C.conj()*C.T*BO*BQ)
    offdiag=~np.eye(model.n,dtype=bool)
    negative=int(np.count_nonzero((norm<0)&offdiag));min_raw=float(norm[offdiag].min())
    norm=np.maximum(norm,0.);np.fill_diagonal(norm,-np.inf)
    partners=np.argmax(norm,axis=1)
    edges=sorted({tuple(sorted((i,int(j))))for i,j in enumerate(partners)})
    nearest=set(fixed_edges(model.shape));overlap=sum(e in nearest for e in edges)
    coords=np.asarray(public['coords']);dist=np.array([np.linalg.norm(coords[i]-coords[j])for i,j in edges])
    man=np.array([abs(i//24-j//24)+abs(i%24-j%24)for i,j in edges])
    pair_array=np.asarray(edges,int);strength=np.array([norm[i,j]for i,j in edges])
    np.savez_compressed(directory/'INTERACTION_GRAPH.npz',partners=partners,pairs=pair_array,
        squared_mixed_response_scores=norm,selected_pair_strength=strength,
        selected_physical_distance=dist,selected_grid_manhattan=man)
    cost['all_dressed_interaction_pair_scores']+=model.n*(model.n-1)//2
    cost['one_partner_argmax_rows']+=model.n
    metadata=dict(selection='percell j!=i argmax ||Oi*Cij*Qj+Oj*Cji*Qi||F²; firstindexbreaksties',
        formula='|Cij|² BOii BQjj+|Cji|² BOjj BQii+2Re(conj(Cij)Cji BOij BQij)',
        selected_unordered_pairs=len(edges),nearest_overlap=overlap,nonlocal_pairs=len(edges)-overlap,
        all_in_previous_nearest_dictionary=overlap==len(edges),
        physical_distance_min=float(dist.min()),physical_distance_median=float(np.median(dist)),physical_distance_max=float(dist.max()),
        grid_manhattan_min=int(man.min()),grid_manhattan_median=float(np.median(man)),grid_manhattan_max=int(man.max()),
        raw_negative_roundoff_scores=negative,raw_min_score=min_raw,score_units='squaredFrobenius mixedsecondresponse, unnormalizedallsourcechannels',
        norm_nonnegative_handling='max(score,0) after exactGramformula; no selectionthreshold',
        zero_score_partner_rows=int(np.count_nonzero(norm[np.arange(model.n),partners]==0)),
        new_resolvent='fresh paid physicalLU+576RHS; no cachedC reuse claim',
        wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu)
    return edges,metadata

if __name__=='__main__':
    with threadpool_limits(limits=1):run('dressed_neighbor_rays',select)
