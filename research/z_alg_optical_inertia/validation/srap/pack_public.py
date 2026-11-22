"""Package only the unmodified released SRAP observations and sampled radii."""

from flint import ctx
import numpy as np

from online import HERE, load_public


if __name__ == '__main__':
    ctx.prec=112
    rows=load_public()
    np.savez_compressed(HERE/'public_observations.npz',
                         ids=np.array([meta['id'] for _,meta in rows]),
                         Y=np.stack([y for y,_ in rows]),
                         original_delta=np.array([meta['stored_delta'] for _,meta in rows]))
