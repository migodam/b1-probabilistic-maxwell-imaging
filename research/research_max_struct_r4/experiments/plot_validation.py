"""Render stored, post-seal endpoints only; no physics, scoring, or selection.

Confirmation truth is accessed only after CONFIRMATION_SCORES.json and its
all18 seal exist. All six independent objects and every frozen primary arm are
shown, including missing/failed endpoints. No development figure is overwritten.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures'
METHODS = ('onehop', 'global', 'native_extended')

def read(path):
    return json.loads(Path(path).read_text())

def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path

def material(path):
    with np.load(resolve(path), allow_pickle=False) as data:
        return data['m'].copy()

def show(ax, values, title):
    n = int(round(np.sqrt(values.size)))
    im = ax.imshow(values.reshape(n, n), origin='lower', vmin=0, vmax=20,
                   cmap='viridis', interpolation='nearest')
    ax.set_title(title, fontsize=10, pad=6)
    ax.set_xticks([])
    ax.set_yticks([])
    return im

def metrics(row):
    return 'NRMSE %.3f | IoU %.3f' % (row['amplitude_normalized_nrmse'], row['otsu_iou'])

def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / (name + '.' + ext), dpi=180)
    plt.close(fig)

def development():
    base = read(ROOT / 'results/DEVELOPMENT_SCORES.json')
    plus = next(r for r in base['rows'] if r['arm'] == 'onehop_absorption')
    phase = read(ROOT / 'results/DEVELOPMENT_PHASE_PI_SCORES_V1.json')
    grid = read(ROOT / 'results/DEVELOPMENT_GRID32_SCORES_V1.json')
    for doc in (phase, grid):
        if not resolve(doc['seal']).is_file():
            raise RuntimeError('Missing development seal')
    pi = phase['rows'][0]
    gridrows = {r['method']: r for r in grid['rows']}
    truth24 = material(ROOT.parent / 'research_max_struct_r2/sealed_truth/obs_b20d5a55677d4eca.npz')
    truth32 = material(ROOT / 'grid32_diagnostic/development_truth/grid32_curved_bar_truth.npz')
    specs = [(material(plus['prediction']), 'Phase + / onehop (24 x 24)\n' + metrics(plus)),
             (material(pi['prediction']), r'Phase $\pi$ / onehop (24 x 24)' + '\n' + metrics(pi)),
             (truth24, 'GT (24 x 24)\nDISCLOSED_DEVELOPMENT'),
             (truth32, 'GT (32 x 32)\nGRID_DIAGNOSTIC'),
             (material(gridrows['global']['prediction']), 'Global (32 x 32)\n' + metrics(gridrows['global'])),
             (material(gridrows['onehop']['prediction']), 'Onehop (32 x 32)\n' + metrics(gridrows['onehop']))]
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 8), layout='constrained')
    for ax, (values, title) in zip(axes.flat, specs):
        im = show(ax, values, title)
    fig.colorbar(im, ax=axes, label='Physical real contrast', shrink=.85, ticks=[0, 5, 10, 15, 20])
    fig.suptitle('DISCLOSED_DEVELOPMENT / GRID_DIAGNOSTIC\nPost-seal stored endpoints; shared contrast [0, 20]', fontsize=13)
    name = 'validation_development_phase_grid'
    save(fig, name)
    return name

def confirmation():
    scorepath = ROOT / 'results/CONFIRMATION_SCORES.json'
    if not scorepath.is_file():
        return None
    scores = read(scorepath)
    sealpath = resolve(scores['seal'])
    if scores.get('status') != 'CONFIRMATION_POST_SEAL_SCORED' or not sealpath.is_file():
        raise RuntimeError('Confirmation scored status and seal required before any truth access')
    sealed = read(sealpath)
    if sealed.get('status') != 'ALL18_DISPOSITIONS_SEALED' or len(sealed['entries']) != 18:
        raise RuntimeError('Complete all18 seal required')
    rows = scores['rows']
    ids = sorted({r['observation_id'] for r in rows})
    if len(ids) != 6 or len(rows) != 18:
        raise RuntimeError('Must retain all6 objects and all18 fixed-arm rows')
    index = {(r['observation_id'], r['method']): r for r in rows}
    endpoints = {(r['observation_id'], r['method']): r for r in sealed['entries']}
    if len(index) != 18 or set(index) != {(i, m) for i in ids for m in METHODS}:
        raise RuntimeError('Incomplete or duplicate primary rows')
    # First confirmation truth access in this script, after all checks above.
    fig, axes = plt.subplots(6, 4, figsize=(13.2, 19.2), layout='constrained')
    for k, ident in enumerate(ids):
        family = index[ident, 'onehop']['family']
        truth = material(ROOT / 'sealed_truth/confirmation' / (ident + '.npz'))
        im = show(axes[k, 0], truth, family + ' / GT\n' + ident)
        for j, method in enumerate(METHODS, 1):
            row = index[ident, method]
            title = family + ' / ' + method + '\n'
            if row.get('gate_evaluation') == 'EVALUATED':
                im = show(axes[k, j], material(endpoints[ident, method]['prediction']), title + metrics(row))
            else:
                axes[k, j].set_facecolor('#eeeeee')
                axes[k, j].text(.5, .5, row['disposition'] + '\nNo scored endpoint', ha='center', va='center', transform=axes[k, j].transAxes)
                axes[k, j].set_title(title + 'NRMSE unavailable | IoU unavailable', fontsize=10)
                axes[k, j].set_xticks([])
                axes[k, j].set_yticks([])
    fig.colorbar(im, ax=axes, label='Physical real contrast', shrink=.7, ticks=[0, 5, 10, 15, 20], pad=.02)
    fig.suptitle('FRESH CONFIRMATION: 6 independent objects / 3 fixed algorithms\nPost-seal material endpoints; all objects and failures retained; shared contrast [0, 20]', fontsize=14)
    name = 'confirmation_six_materials_postseal'
    save(fig, name)
    return name

def main():
    tick = time.perf_counter()
    OUT.mkdir(exist_ok=True)
    names = [development()]
    fresh = confirmation()
    if fresh:
        names.append(fresh)
    receipt = dict(status='RENDERED', rendered_utc=datetime.now(timezone.utc).isoformat(),
        truth_access='POST_SEAL_ONLY', confirmation_truth_access=fresh is not None,
        confirmation_status='RENDERED_ALL6_ALL18' if fresh else 'NOT_RUN_SCORES_NOT_AVAILABLE',
        figures=[str(OUT / (n + '.' + ext)) for n in names for ext in ('png', 'pdf')],
        new_physical_LU=0, new_complex_RHS=0, new_metrics=0, candidate_selection=False,
        existing_development_materials_preserved=True, wall_seconds=time.perf_counter()-tick)
    (ROOT / 'audit/VALIDATION_PLOT_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    main()
