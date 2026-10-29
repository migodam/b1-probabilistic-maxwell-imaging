"""Ordinary public binary 5x7 digit strokes in 7x7 canvases.

Each digit is a separate static material-support world, never a class label to
be inferred. Rows run TOP TO BOTTOM; columns run LEFT TO RIGHT. Physical cell
size, position, frequencies and inversion grid are deliberately not assigned.
"""
from pathlib import Path
import numpy as np

PATTERN_ROWS = {
    1: ('00100','01100','00100','00100','00100','00100','01110'),
    2: ('01110','10001','00001','00010','00100','01000','11111'),
    3: ('11110','00001','00001','01110','00001','00001','11110'),
    4: ('00010','00110','01010','10010','11111','00010','00010'),
    5: ('11111','10000','10000','11110','00001','00001','11110'),
    6: ('01110','10000','10000','11110','10001','10001','01110'),
    7: ('11111','00001','00010','00100','01000','01000','01000'),
    8: ('11111','10001','10001','11111','10001','10001','11111'),
    9: ('01110','10001','10001','01111','00001','00001','01110'),
}
DIGIT_ORDER = (8, 1, 2, 3, 4, 5, 6, 7, 9)
DISPLAY_ORDER = tuple(range(1, 10))


def get_pattern(digit, padded=True):
    """Return a fresh binary uint8 array: (7,7) padded, or (7,5) stroke."""
    if digit not in PATTERN_ROWS:
        raise ValueError('Only public digits 1..9 are defined')
    stroke = np.asarray([[int(c) for c in row] for row in PATTERN_ROWS[digit]], dtype=np.uint8)
    if not padded:
        return stroke
    canvas = np.zeros((7, 7), dtype=np.uint8)
    canvas[:, 1:6] = stroke
    return canvas


def digit_mask(digit, padded=True):
    return get_pattern(digit, padded=padded)


def export():
    import json, time
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cpu, wall = time.process_time(), time.perf_counter()
    out = Path(__file__).resolve().parent
    arrays = {}
    for digit in DISPLAY_ORDER:
        stroke, canvas = get_pattern(digit, False), get_pattern(digit)
        arrays[f'digit_{digit}'] = canvas
        np.savez_compressed(out / f'digit_{digit}.npz', stroke_5x7=stroke, canvas_7x7=canvas, mask=canvas)
        fig, ax = plt.subplots(figsize=(2.7, 2.9))
        ax.imshow(canvas, cmap='gray_r', vmin=0, vmax=1, interpolation='nearest', origin='upper')
        ax.set_title(f'Digit {digit} — independent support', fontsize=10)
        ax.set_xticks(np.arange(-.5, 7, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 7, 1), minor=True)
        ax.grid(which='minor', color='#cccccc', linewidth=.55)
        ax.tick_params(which='both', left=False, bottom=False, labelleft=False, labelbottom=False)
        fig.tight_layout()
        fig.savefig(out / f'digit_{digit}.png', dpi=160, bbox_inches='tight')
        plt.close(fig)
    np.savez_compressed(out / 'DIGIT_PATTERNS_1_TO_9.npz', **arrays)
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 7.8))
    for digit, ax in zip(DISPLAY_ORDER, axes.ravel()):
        ax.imshow(arrays[f'digit_{digit}'], cmap='gray_r', vmin=0, vmax=1, interpolation='nearest', origin='upper')
        ax.set_title(f'{digit}' + ('  (first exploration)' if digit == 8 else ''), fontsize=11)
        ax.set_xticks(np.arange(-.5, 7, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 7, 1), minor=True)
        ax.grid(which='minor', color='#cccccc', linewidth=.45)
        ax.tick_params(which='both', left=False, bottom=False, labelleft=False, labelbottom=False)
    fig.suptitle('Public digit supports: nine separate static worlds', fontsize=13)
    fig.text(.5, .02, '5×7 binary strokes; 7×7 canvas with blank left/right columns. Physical scale unset.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .96))
    fig.savefig(out / 'DIGITS_1_TO_9.png', dpi=180, bbox_inches='tight')
    plt.close(fig)
    metadata = dict(
        schema='PUBLIC_DIGIT_SUPPORTS_5X7_IN_7X7_V1',
        world_semantics='Each digit is one independent static material-support geometry; not a classification task or label-selected template dictionary for inversion.',
        construction='Explicit ordinary binary stroke patterns in this module; no private inputs or external font assets.',
        row_order='top_to_bottom', column_order='left_to_right', dtype='uint8', foreground=1, background=0,
        stroke_shape=[7,5], canvas_shape=[7,7], blank_columns=[0,6],
        display_order=list(DISPLAY_ORDER), first_exploration_digit=8, exploration_order=list(DIGIT_ORDER),
        physical_scale=None, physical_position=None, GHz=None, inversion_grid=None,
        physical_configuration_owner='Eval/root freezes all physical scales and grids separately',
        patterns={str(d):dict(rows_5x7=list(PATTERN_ROWS[d]), foreground_cells=int(arrays[f'digit_{d}'].sum()),
                              npz=f'digit_{d}.npz', png=f'digit_{d}.png', physics_status='NOT_RUN_PHYSICS') for d in DISPLAY_ORDER},
        note='Digit8 contains two enclosed background holes in this binary topology; physical resolution remains untested.',
        generation_receipt=dict(full_fe=0,rhs=0,process_cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
    (out / 'PATTERN_METADATA.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(out),first_exploration_digit=8,physical_cost=dict(full_fe=0,rhs=0),files='nine raw NPZ + nine PNG, aggregate NPZ, contact PNG, metadata')))


if __name__ == '__main__':
    export()
