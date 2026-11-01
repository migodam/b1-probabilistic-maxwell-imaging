"""Plot the existing two-cell algebraic example; do not alter its JSON evidence."""
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/fbg_r3_mpl_cache')
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HERE = Path(__file__).resolve().parent
case = json.loads((HERE / 'branches_checks.json').read_text())['isolated_roots']
gparts = np.array(case['G'])
G = gparts[:, :, 0] + 1j * gparts[:, :, 1]
E, S = np.array(case['E']), np.array(case['S'])
Y = complex(*case['y'])
ma = np.array(case['roots'][1]['m'])
mb = np.array(case['roots'][0]['m'])

def current(m):
    return np.linalg.solve(np.eye(2) - m[:, None] * G, m * E)

ja, jb = current(ma), current(mb)
qa, qb = E + G @ ja, E + G @ jb
tgrid = np.linspace(0, 1, 401)
state_linear, state_profile = [], []
data_linear, data_profile, data_chord = [], [], []
profile_materials, physical_min_sv, chord_check, profile_check = [], [], [], []
for t in tgrid:
    jt = (1 - t) * ja + t * jb
    qt = E + G @ jt
    mt = (1 - t) * ma + t * mb
    mp = np.real(qt.conj() * jt) / np.abs(qt) ** 2
    dl = jt - mt * qt
    dp = jt - mp * qt
    state_linear.append(np.linalg.norm(dl))
    state_profile.append(np.linalg.norm(dp))
    data_linear.append(abs(S @ current(mt) - Y) / abs(Y))
    data_profile.append(abs(S @ current(mp) - Y) / abs(Y))
    data_chord.append(abs(S @ jt - Y) / abs(Y))
    profile_materials.append(mp)
    physical_min_sv.extend([
        np.linalg.svd(np.eye(2) - m[:, None] * G, compute_uv=False)[-1]
        for m in (mt, mp)
    ])
    chord_check.append(np.linalg.norm(dl - t * (1 - t) * (mb - ma) * (qb - qa)))
    numerator = (np.abs(qa) ** 2 * np.abs(qb) ** 2
                 - np.real(qa.conj() * qb) ** 2)
    angle = t ** 2 * (1 - t) ** 2 * (mb - ma) ** 2 * numerator / np.abs(qt) ** 2
    profile_check.append(np.linalg.norm(np.abs(dp) ** 2 - angle))

state_linear, state_profile, data_linear, data_profile, data_chord = map(
    np.asarray, (state_linear, state_profile, data_linear, data_profile, data_chord)
)
profile_materials = np.array(profile_materials)
blue, orange, ink = '#234f83', '#c96b23', '#26313f'
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10.5,
    'axes.labelsize': 11, 'axes.titlesize': 12,
    'axes.edgecolor': '#9ba6b2', 'axes.labelcolor': ink,
    'xtick.color': '#505c69', 'ytick.color': '#505c69',
    'text.color': ink, 'grid.color': '#dfe5eb', 'grid.linewidth': .7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'savefig.facecolor': 'white',
})
fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.7))
fig.subplots_adjust(left=.075, right=.975, bottom=.20, top=.76, wspace=.26)
fig.suptitle('Two-cell algebraic VIE — not an imaging benchmark',
             x=.075, y=.965, ha='left', fontsize=16, fontweight='bold')
fig.text(.075, .907,
         'Same-data endpoints: (2, 7) and (4.54062, 2.02148).  '
         'The current follows their straight chord.', fontsize=11, color='#526171')

left, right = axes
for ax in axes:
    ax.set_xlim(0, 1)
    ax.set_xlabel('Interpolation parameter t')
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.grid(True, axis='both', zorder=0)
    ax.set_axisbelow(True)

left.plot(tgrid, state_linear, color=blue, lw=2.5, label='Linear material', zorder=3)
left.plot(tgrid, state_profile, color=orange, lw=2.5, label='Profiled real material', zorder=4)
left.set_ylim(-.009, .48)
left.set_title('A   State defect of the lifted current', loc='left', pad=13, fontweight='bold')
left.set_ylabel(r'$\|j_t-X_t(E+Gj_t)\|_2$')
left.legend(loc='upper left', frameon=False, fontsize=9.5)
left.text(.51, .433, f'Peak {state_linear.max():.3f}', color=blue,
          ha='center', fontsize=10)

right.plot(tgrid, data_linear, color=blue, lw=2.5, label='Linear material: physical solve', zorder=3)
right.plot(tgrid, data_profile, color=orange, lw=2.5, label='Profiled material: physical solve', zorder=4)
right.plot(tgrid, data_chord, color=ink, lw=1.4, ls=(0, (5, 3)),
           label='Current chord: analytically zero', zorder=5)
right.set_ylim(-.002, .101)
right.set_title('B   Data residual after re-solving Maxwell', loc='left', pad=13, fontweight='bold')
right.set_ylabel(r'$|F(m_t)-Y|\,/\,|Y|$')
right.legend(loc='upper left', frameon=False, fontsize=9.0)
right.text(.79, .0900, f'Peak {data_linear.max():.4f}', color=blue,
           ha='center', fontsize=10)

# The insets retain exact values and expose the smaller orange curves.
for ax, values, upper in [(left, state_profile, .0178), (right, data_profile, .0043)]:
    sub = ax.inset_axes([.29, .18, .51, .31])
    sub.set_facecolor('#fffaf5')
    sub.plot(tgrid, values, color=orange, lw=1.9)
    sub.set_xlim(0, 1)
    sub.set_ylim(0, upper)
    sub.set_xticks([0, .5, 1])
    sub.yaxis.set_major_locator(MaxNLocator(3))
    sub.tick_params(labelsize=8, length=2.5)
    sub.grid(True, color='#eadfD3', lw=.55)
    sub.set_title(f'Profile detail: peak {values.max():.4f}',
                  fontsize=8.5, color=orange, pad=4)
    for spine in sub.spines.values():
        spine.set_color('#d8bda3')

fig.text(.075, .095,
         r'The chord has $Sj_t=Y$ throughout.  '
         'Re-solving Maxwell for its material readouts reveals nonzero data residuals between the roots.',
         fontsize=10, color='#465565')
fig.text(.075, .044,
         'Insets show the profiled curves on their own scales.  '
         'Both endpoints are physical roots; no new reconstruction experiment is shown.',
         fontsize=9.5, color='#667383')
fig.savefig(HERE / 'branch_corridor.png', dpi=240)
plt.close(fig)

print(json.dumps({
    'figure': str(HERE / 'branch_corridor.png'),
    'max_state_linear': float(state_linear.max()),
    'max_state_profile': float(state_profile.max()),
    'max_relative_physical_data_linear': float(data_linear.max()),
    'max_relative_physical_data_profile': float(data_profile.max()),
    'max_relative_chord_data': float(data_chord.max()),
    'minimum_resolvent_singular_value_on_both_material_paths': float(min(physical_min_sv)),
    'minimum_real_field_cross_product': float(np.real(qa.conj()*qb).min()),
    'profile_material_in_endpoint_box': bool(np.all(profile_materials >= np.minimum(ma, mb)-1e-12)
                                             and np.all(profile_materials <= np.maximum(ma, mb)+1e-12)),
    'max_chord_identity_absolute_error': float(max(chord_check)),
    'max_profile_angle_squared_defect_absolute_error': float(max(profile_check)),
    'source_json_modified': False,
}, indent=2))
