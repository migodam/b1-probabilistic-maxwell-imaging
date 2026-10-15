"""One continuous phase-field path, software representation diagnostic only."""
from pathlib import Path
import json
import numpy as np
from scipy.special import expit
from g_archive import structure_features, topology_summary, material_distance
from g_budget import GBudget, dump

ROOT = Path(__file__).resolve().parents[1]


def run():
    with GBudget(ROOT).job('g-e0-continuous-representation-v1', stage='E0'):
        points = (np.arange(24) + .5) / 24 - .5
        y, x = np.meshgrid(points, points, indexing='ij')
        r2 = x*x + y*y
        parameters = np.linspace(-.02, .07, 19)
        fields = []
        for z in parameters:
            # Location and outer disk are declared toy constants, not a target
            # truth or a fitted template. No Maxwell data enter this example.
            amplitude = 1.2 * expit((.4**2-r2)/.002) * expit((r2-z)/.002)
            fields.append(amplitude.astype(np.complex128))
        fields = np.array(fields)
        diagnostics = [topology_summary(structure_features(f, np.zeros((24, 24)))) for f in fields]
        dump(ROOT/'results/CONTINUOUS_REPRESENTATION.json', dict(
            status='SOFTWARE_REPRESENTATION_ONLY_NOT_MAXWELL_FEASIBLE_DISCOVERY',
            formula='1.2 sigmoid((.4^2-r^2)/.002) sigmoid((r^2-z)/.002)',
            parameters=parameters, topology=[d['holes_ge6'] for d in diagnostics],
            adjacent_material_distances=[material_distance(a, b) for a, b in zip(fields[:-1], fields[1:])],
            supplied_hole_count_labels=False, true_material_read=False,
            limitations='A fixed centered toy path illustrates topology changes under continuous parameters; it does not represent arbitrary shapes, posterior modes or Maxwell-feasible worlds',
            new_states=0, new_RHS=0))
        np.savez_compressed(ROOT/'results/CONTINUOUS_REPRESENTATION.npz', z=parameters, fields=fields)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        figure, axes = plt.subplots(1, 5, figsize=(12, 3), layout='constrained')
        for ax, i in zip(axes, [0, 4, 8, 12, 18]):
            ax.imshow(fields[i].real, origin='lower', vmin=0, vmax=1.2, cmap='viridis')
            ax.set_title(f'z={parameters[i]:.3f}; holes={diagnostics[i]["holes_ge6"]}', fontsize=9)
            ax.set_axis_off()
        figure.suptitle('Continuous parameter, changing support: software-only; no Maxwell validation', fontsize=10)
        (ROOT/'figures').mkdir(exist_ok=True)
        figure.savefig(ROOT/'figures/continuous_representation.png', dpi=130)
        plt.close(figure)
        source = json.loads((ROOT/'results/E0_TRAJECTORY_AUDIT.json').read_text())
        methods = list(source['radial']['summary'])
        fig, ax = plt.subplots(figsize=(9, 4), layout='constrained')
        for j, key in enumerate(['raw_proposals','corrected_endpoints','full_valid_endpoints','historical_archive']):
            values = [source['radial']['summary'][m]['mean_counts_by_stage'][key]['original_material_distinct'] for m in methods]
            ax.bar(np.arange(len(methods))+(j-1.5)*.19, values, width=.19, label=key)
        ax.set_xticks(np.arange(len(methods)), methods)
        ax.set_ylabel('Material representatives per exposed parent')
        ax.set_title('Saved radial trajectories: proposal diversity is not endpoint diversity')
        ax.legend(fontsize=8)
        fig.savefig(ROOT/'figures/E0_stage_attribution.png', dpi=130)
        plt.close(fig)


if __name__ == '__main__':
    run()
