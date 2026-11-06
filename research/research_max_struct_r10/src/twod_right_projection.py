"""Root-prescribed exposed diagnostic: keep G P output, no left projection.

The existing real-material denominator solver/readout is reused unchanged.
No GT, scoring, source-current change, new loss, beta or regularization.
"""
from pathlib import Path
import json,time
import numpy as np
from research.research_max_struct_r10.src.lowrank_denominator import (
    prepare_bases,LowRankDenominator,matrix_material_readout)


def prepare_green_right(model,config):
    started=time.process_time()
    bases,info=prepare_bases(model.G,model.S,.1)
    operator=LowRankDenominator(model.G,bases['green_strong'],projection='right')
    supplied=dict(config,matrix_denominator=operator,matrix_inner_maxiter=8)
    out=Path(config['receipt_path']).parent
    info.update(basis='green_strong',projection='right',matrix_inner_maxiter=8,
        actual_denominator_rank=operator.rank,
        setup_complete_cpu_seconds=time.process_time()-started,
        status='EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC',truth_access='NONE',
        change='K=G U U* instead of U U* G U U*; all other solver settings inherited')
    (out/'DENOMINATOR_SETUP.json').write_text(json.dumps(info,indent=2)+'\n')
    np.savez_compressed(out/'DENOMINATOR_BASIS.npz',U=operator.U,B=operator.B)
    return supplied
