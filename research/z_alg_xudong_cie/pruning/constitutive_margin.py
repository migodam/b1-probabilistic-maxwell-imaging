"""How far a supported data-fitting current must miss its best shared physical material."""
import json
import numpy as np
from flint import arb, ctx
from scipy.fft import fft2, ifft2
from threadpoolctl import threadpool_limits

from native_box import HERE, ball, clock, elapsed, frobenius_squared, hi, lo
from site_driver import DEFAULT_PUBLIC


def apply_green(current, mask, symbol):
    field = np.zeros((96, 96, current.shape[1]), dtype=complex)
    field.reshape(-1, current.shape[1])[mask] = current
    padded = fft2(field, s=symbol.shape, axes=(0, 1))
    image = ifft2(symbol[:, :, None] * padded, axes=(0, 1))[:96, :96]
    return image.reshape(-1, current.shape[1])[mask]


def main():
    start = clock()
    ctx.prec = 128
    online = json.loads((HERE / 'site_online.json').read_text())['predictions']
    reciprocal = json.loads((HERE / 'reciprocal_control.json').read_text())['rows']
    fits = {(r['observation'], tuple(r['tile'])): r for r in reciprocal
            if r['strict_fit'] and r['original_decision'] == 'PRUNE'}
    rows, sensors = [], {}
    grid = np.indices((96, 96)) // 24
    with np.load(DEFAULT_PUBLIC/'geometry96.npz') as geometry, np.load(HERE/'reciprocal_witnesses.npz') as z:
        symbol = geometry['spectrum']
        g = max((arb(float(x.real))**2 + arb(float(x.imag))**2).upper()
                for x in symbol.ravel()).sqrt().upper()
        for name, view in online.items():
            aperture = view['aperture']
            B, S = geometry[aperture + '__B'], geometry[aperture + '__S']
            for entry in view['tile_scores']:
                primary, tile = entry['primary'], tuple(entry['tile'])
                if primary['decision'] != 'PRUNE':
                    continue
                tile_proof = primary.get('translated_tile_check', primary)
                C = arb(tile_proof['current_norm_factor_upper'])
                mask = ((grid[0] == tile[0]) & (grid[1] == tile[1])).ravel()
                if (aperture, tile) not in sensors:
                    sensors[aperture, tile] = frobenius_squared(ball(S[:, mask])).sqrt().upper()
                s = sensors[aperture, tile]
                margin = arb(primary['residual_lower']) - arb(view['radius'])
                lower = lo((margin / (s * (1 + C * g))).lower())
                row = dict(observation=name,tile=tile,data_gap_lower=lo(margin),
                           universal_state_defect_lower=lower,receiver_norm_upper=hi(s))
                record = fits.get((name, tile))
                if record is not None:
                    q, x = z[record['witness_basis']], z[record['witness_matrix']]
                    j = record['current_gain_C'] * q @ x @ (q.T @ B[mask])
                    e = B[mask] + apply_green(j, mask, symbol)
                    power = np.sum(abs(e)**2, axis=1)
                    cross = np.sum(e.conj() * j, axis=1).real
                    m = np.clip(np.divide(cross, power, out=np.zeros_like(cross), where=power>0), 0, 4)
                    row['observed_relative_defect_numerical'] = float(np.linalg.norm(j-m[:, None]*e)/np.linalg.norm(j))
                    row['observed_absolute_defect_numerical'] = float(np.linalg.norm(j-m[:, None]*e))
                rows.append(row)
    if len(rows) != 240 or sum('observed_absolute_defect_numerical' in x for x in rows) != len(fits):
        raise ValueError('Stored support or transfer certificates do not cover the stated cohort.')
    data = dict(scope='Uniform strict lower for any data-compatible current on each ruled-out support. Separate point diagnostics on saved reciprocal fits use binary64.',
                native_green_norm_upper=hi(g), rows=rows, **elapsed(start))
    (HERE / 'constitutive_margin.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    for key in ('universal_state_defect_lower','observed_relative_defect_numerical'):
        values=[r[key] for r in rows if key in r]
        print(key, len(values), np.quantile(values,[0,.5,1]))
    print(elapsed(start))


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
