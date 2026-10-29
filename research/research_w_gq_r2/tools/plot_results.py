"""Plot saved R2 evidence only. Never imports a solver or runs a PDE.

Run through tools/run_job.py so plotting CPU/RSS is charged once.  Raw plot
tables, source paths and limitations accompany every figure in figures/.
"""
from pathlib import Path
import argparse
import csv
import json
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.patches import Ellipse

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT/'figures'
DATA = FIG/'data'
DELEGATED = ROOT/'research/delegated/plots'
COLORS = ['#246a9c', '#d77a24', '#319477', '#8a5da2']
REGISTRY = {}


def read_json(path):
    return json.loads((ROOT/path).read_text())


def read_csv(path):
    with (ROOT/path).open() as stream:
        return list(csv.DictReader(stream))


def raw_csv(name, rows):
    rows = list(rows)
    if not rows:
        raise ValueError('empty raw table: '+name)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    path = DATA/(name+'.csv')
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return str(path.relative_to(ROOT))


def save(fig, name, title, sources, raw, explanation, scope):
    fig.suptitle(title, fontsize=14, fontweight='semibold')
    fig.savefig(FIG/(name+'.png'), dpi=220, facecolor='white', bbox_inches='tight')
    fig.savefig(FIG/(name+'.pdf'), facecolor='white', bbox_inches='tight')
    plt.close(fig)
    REGISTRY[name] = dict(status='PLOTTED_SAVED_DATA', png='figures/'+name+'.png',
        pdf='figures/'+name+'.pdf', raw_csv=raw, source_files=sources,
        explanation_zh=explanation, scope_zh=scope)


def decorate(ax, xlabel=None, ylabel=None, log=False):
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale('log')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=.18, linewidth=.6)
    ax.set_axisbelow(True)


def micro_materials(rows):
    name = 'F01_MICRO_MATERIAL_MAPS'
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6), layout='constrained')
    values = []
    sources = []
    for p in range(2):
        walls = []
        for tex in range(2):
            rel = f'results/micro/WALL_p{p}_n16_t{tex}.npz'
            sources.append(rel)
            with np.load(ROOT/rel) as saved:
                walls.append({k: saved[k].copy() for k in ('x', 'y', 'eps', 'ell')})
        x, y = walls[0]['x'], walls[0]['y']
        for tex, wall in enumerate(walls):
            axes[p, tex].pcolormesh(x, y, wall['eps'], cmap='viridis',
                vmin=1., vmax=4.8, shading='flat', rasterized=True)
            axes[p, tex].set_title(f'Parent {p}, texture {"AB"[tex]}')
            for iy, yc in enumerate((y[:-1]+y[1:])/2):
                for ix, xc in enumerate((x[:-1]+x[1:])/2):
                    values.append(dict(parent=p, texture=tex, n=16, eta=1/16,
                        ix=ix, iy=iy, x=xc, y=yc, eps_real=wall['eps'][iy, ix],
                        loss_parameter=wall['ell'][iy, ix], source= sources[-2+tex]))
        delta = walls[0]['eps']-walls[1]['eps']
        difference = axes[p, 2].pcolormesh(x, y, delta, cmap='RdBu_r',
            vmin=-1.6, vmax=1.6, shading='flat', rasterized=True)
        axes[p, 2].set_title(f'Parent {p}: A minus B')
        for ax in axes[p]:
            ax.set_xlim(-.85, .85)
            ax.set_ylim(.48, 1.27)
            ax.set_aspect('equal')
            ax.set_xlabel(r'$x/\lambda_0$')
            ax.set_ylabel(r'$y/\lambda_0$')
    fig.colorbar(axes[0, 0].collections[0], ax=axes[:, :2], shrink=.8,
        label=r'Relative permittivity $\epsilon_\infty$')
    fig.colorbar(difference, ax=axes[:, 2], shrink=.8,
        label=r'$\epsilon_{\infty,A}-\epsilon_{\infty,B}$')
    table = raw_csv(name, values)
    save(fig, name, 'Two resolved microstructure walls at eta / lambda0 = 1/16',
         sources, [table], '两个宏观parent的实际材料A/B及其差，采用保存的墙单元边界和材料值。',
         '颜色为真实离散材料，不表示反演恢复；所有父世界及两纹理均显示，空气为1。')


def micro_scale(rows):
    name = 'F02_MICRO_SCALE_AND_ALL_PROBES'
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
    summary = []
    for parent in range(2):
        points = []
        for n in (16, 12, 8):
            rr = [r for r in rows if r['parent'] == parent and r['n'] == n]
            common = [r for r in rr if r['f'] == 1 and r['target'] < 4]
            held = [r for r in rr if r['f'] == 1 and r['role'] == 'heldout']
            row = dict(parent=parent, n=n, eta=1/n,
                material_RMS_occupied=rr[0]['material_RMS_occupied'],
                material_RMS_public=rr[0]['material_RMS_public'],
                common_four_A0_max=max(r['d_A0'] for r in common),
                common_four_A1_max=max(r['d_A1'] for r in common),
                all_available_A1_max=max(r['d_A1'] for r in rr),
                all_available_count=len(rr),
                heldout_f1_count=len(held),
                heldout_f1_A1_max=max(r['d_A1'] for r in held))
            points.append(row)
            summary.append(row)
        x = [r['eta'] for r in points]
        color = COLORS[parent]
        axes[0].plot(x, [r['material_RMS_occupied'] for r in points],
            'o-', color=color, label=f'P{parent}, occupied')
        axes[0].plot(x, [r['material_RMS_public'] for r in points],
            'o--', color=color, alpha=.7, label=f'P{parent}, public ROI')
        axes[1].plot(x, [r['common_four_A1_max'] for r in points],
            'o-', color=color, label=f'P{parent}: matched 4, f0')
        axes[1].plot(x, [r['all_available_A1_max'] for r in points],
            'D:', color=color, alpha=.8, label=f'P{parent}: all saved probes')
        rr = sorted([r for r in rows if r['parent'] == parent and r['n'] == 16
                     and r['f'] == 1], key=lambda r:r['target'])
        axes[2].plot([r['target'] for r in rr], [r['d_A1'] for r in rr],
            'o-', color=color, ms=4, label=f'Parent {parent}')
    axes[0].axhline(.5, color='#777', ls=':', label='Prespecified material scale 0.5')
    for ax in axes[1:]:
        ax.axhline(1., color='#555', ls='--', lw=1, label='Pair-distance scale 1')
    for ax in axes[:2]:
        ax.set_xticks([1/16, 1/12, 1/8], ['1/16', '1/12', '1/8'])
    axes[2].axvspan(-.3, 1.5, color='#aaa', alpha=.13)
    axes[2].text(.04, .89, 'Shaded: train targets 0, 1', transform=axes[2].transAxes,
        va='top', fontsize=8)
    axes[0].set_title('Material contrast remains visible')
    axes[1].set_title('Four-source response: actual maxima')
    axes[2].set_title('Finest scale: all 14 probes at f0')
    decorate(axes[0], r'$\eta/\lambda_0$', 'Material RMS difference')
    decorate(axes[1], r'$\eta/\lambda_0$', 'Whitened pair distance', log=True)
    decorate(axes[2], 'Frozen target index', 'Whitened pair distance', log=True)
    for ax in axes:
        ax.legend(fontsize=7.5, loc='best')
    tables = [raw_csv(name+'_RAW', rows), raw_csv(name+'_SUMMARY', summary)]
    save(fig, name, 'Material difference and finite-task full-wave response difference',
        ['results/micro/RESULTS.json'], tables,
        '左图显示真实材料RMS；中图分别显示所有尺度共同四探针最大值及所有已保存探针最大值；右图完整列出最细尺度14个f0目标。',
        '各尺度全部探针数不同：n=8/12只运行4个f0探针及两个单目标频率压力点，n=16有14个f0探针。全保存最大值不是整个连续目标族上界。')


def micro_validation():
    name = 'F03_INDEPENDENT_GRID_COMPARISON'
    rows = read_json('results/micro_validation/RESULTS.json')
    paired = []
    endpoints = []
    for parent in range(2):
        for f in (1., 1.15):
            rr = [r for r in rows if r['parent'] == parent and r['f'] == f]
            if len(rr) != 2:
                continue
            rr.sort(key=lambda r:r['texture'])
            saved = []
            for r in rr:
                rel = f'results/micro_validation/p{parent}_f{f}_t{r["texture"]}.npz'
                with np.load(ROOT/rel) as z:
                    saved.append((z['coarse'].copy(), z['reference'].copy()))
                endpoints.append(dict(parent=parent, f=f, target=r['target'],
                    texture=r['texture'], endpoint_error_white=r['endpoint_error_white'],
                    reference_ndofs=r['ndofs'], source=rel))
            coarse = np.linalg.norm(saved[0][0]-saved[1][0])/.01
            reference = np.linalg.norm(saved[0][1]-saved[1][1])/.01
            pair_change = np.linalg.norm((saved[0][1]-saved[1][1])-
                                        (saved[0][0]-saved[1][0]))/.01
            paired.append(dict(parent=parent, f=f, target=rr[0]['target'],
                coarse_pair_distance=coarse, reference_pair_distance=reference,
                pair_difference_discrepancy=pair_change,
                empirical_endpoint_error_sum=sum(r['endpoint_error_white'] for r in rr)))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), layout='constrained')
    x = np.arange(len(paired))
    labels = [f'P{r["parent"]}, f={r["f"]:g}\nprobe {r["target"]}' for r in paired]
    axes[0].bar(x-.17, [r['coarse_pair_distance'] for r in paired], width=.34,
        color=COLORS[0], label='Primary grid')
    axes[0].bar(x+.17, [r['reference_pair_distance'] for r in paired], width=.34,
        color=COLORS[1], label='Independent refinement')
    axes[0].set_xticks(x, labels)
    axes[0].legend(fontsize=8)
    for tex in range(2):
        axes[1].bar(x+(tex-.5)*.34,
            [next(e['endpoint_error_white'] for e in endpoints if e['parent']==r['parent']
                  and e['f']==r['f'] and e['texture']==tex) for r in paired],
            width=.34, color=COLORS[tex], label=f'Texture {"AB"[tex]}')
    axes[1].axhline(.1, color='#555', ls='--', label='Qualification scale 0.1')
    axes[1].set_xticks(x, labels)
    axes[1].legend(fontsize=8)
    axes[0].set_title('Same selected probe, two physical grids')
    axes[1].set_title('Endpoint discretization differences')
    decorate(axes[0], ylabel='Four-source whitened pair distance')
    decorate(axes[1], ylabel='Endpoint whitened response difference', log=True)
    tables=[raw_csv(name+'_PAIRS', paired), raw_csv(name+'_ENDPOINTS', endpoints)]
    save(fig, name, 'Independent refinement at each parent / frequency worst probe',
        ['results/micro_validation/RESULTS.json']+[r['source'] for r in endpoints], tables,
        '每个parent及频率在自身argmax处独立细化，并直接比较保存的两端响应及pair差；没有借用其他probe的误差。',
        '经验网格差及端点差和不是连续PDE误差证书，也不覆盖未复核的全部目标。')


def reference_subtriangles(n=3):
    triangles = []
    for i in range(n):
        for j in range(n-i):
            triangles.append(np.array(((i,j),(i+1,j),(i,j+1)), float)/n)
            if i+j < n-1:
                triangles.append(np.array(((i+1,j),(i+1,j+1),(i,j+1)), float)/n)
    return np.asarray(triangles)


def saved_p3_samples(path, source=0):
    """Rebuild only DOF topology, then evaluate saved P3 polynomials locally.

    The physical field coefficients are saved free DOFs, not mesh vertices.
    No stiffness/mass assembly, factorization, solve or global locator is used.
    """
    from skfem import Basis, MeshTri, ElementTriP3
    with np.load(path) as z:
        p, t, u = z['mesh_p'].copy(), z['mesh_t'].copy(), z['u'].copy()
        target = json.loads(z['target'].item())
    mesh = MeshTri(p, t)
    if not np.array_equal(mesh.t, t):
        raise ValueError('saved/reconstructed element order changed')
    element = ElementTriP3()
    basis = Basis(mesh, element, intorder=1)
    boundary = basis.get_dofs().all()
    free = np.setdiff1d(np.arange(basis.N), boundary)
    if u.shape[0] != len(free):
        raise ValueError('saved free P3 DOF count differs from reconstructed topology')
    full = np.zeros(basis.N, complex)
    full[free] = u[:, source]
    centers = p[:, t].mean(axis=1)
    selected = np.flatnonzero((centers[0] >= -1.6) & (centers[0] <= 1.6)
        & (centers[1] >= .35) & (centers[1] <= 2.9))
    subdiv = reference_subtriangles(3)
    points = subdiv.mean(axis=1).T
    phi = np.array([element.lbasis(points, j)[0] for j in range(basis.Nbfun)])
    values = phi.T@full[basis.element_dofs[:, selected]]
    a = p[:, t[0, selected]]
    b = p[:, t[1, selected]]-a
    c = p[:, t[2, selected]]-a
    polygons = np.concatenate([(a[None]+r[:, 0, None, None]*b[None]
        +r[:, 1, None, None]*c[None]).transpose(2,0,1) for r in subdiv])
    centroids = polygons.mean(axis=1)
    return polygons, values.reshape(-1), centroids, target, dict(
        saved_free_dofs=int(u.shape[0]), reconstructed_total_dofs=int(basis.N),
        sampled_elements=len(selected), samples_per_element=len(subdiv),
        source=source, field='scattered field u_sc',
        evaluation='saved P3 coefficients at known-element subtriangle centroids',
        rendering='flat subtriangles; actual P3 polynomial samples; no PDE')


def micro_fields():
    name = 'F04_SAVED_P3_COMPLEX_SCATTERED_FIELD'
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.6), layout='constrained')
    samples = []
    schemas = []
    sources = []
    plotted = []
    for tex in range(2):
        rel = f'results/micro/p0_n16_t{tex}_f1_q1.npz'
        sources.append(rel)
        polygons, values, centers, target, meta = saved_p3_samples(ROOT/rel)
        schemas.append(dict(texture=tex, source_file=rel, **meta))
        plotted.append((polygons, values, centers, target))
        for i, (xy, v) in enumerate(zip(centers, values)):
            samples.append(dict(parent=0, texture=tex, eta=1/16, f=1, source=0,
                sample=i, x=xy[0], y=xy[1], u_sc_real=v.real, u_sc_imag=v.imag))
    realmax = max(float(np.max(abs(q[1].real))) for q in plotted)
    imagmax = max(float(np.max(abs(q[1].imag))) for q in plotted)
    for tex, (polygons, values, centers, target) in enumerate(plotted):
        wallfile = ROOT/f'results/micro/WALL_p0_n16_t{tex}.npz'
        with np.load(wallfile) as w:
            wx, wy, occupied = w['x'].copy(), w['y'].copy(), w['eps']>1
        for component, vmax in enumerate((realmax, imagmax)):
            ax = axes[tex, component]
            val = values.real if component == 0 else values.imag
            collection = PolyCollection(polygons, array=val, cmap='RdBu_r',
                clim=(-vmax, vmax), linewidth=0, antialiased=False, rasterized=True)
            ax.add_collection(collection)
            edges=[]
            for iy,ix in np.argwhere(occupied):
                if ix==0 or not occupied[iy,ix-1]:
                    edges.append(((wx[ix],wy[iy]),(wx[ix],wy[iy+1])))
                if ix==occupied.shape[1]-1 or not occupied[iy,ix+1]:
                    edges.append(((wx[ix+1],wy[iy]),(wx[ix+1],wy[iy+1])))
                if iy==0 or not occupied[iy-1,ix]:
                    edges.append(((wx[ix],wy[iy]),(wx[ix+1],wy[iy])))
                if iy==occupied.shape[0]-1 or not occupied[iy+1,ix]:
                    edges.append(((wx[ix],wy[iy+1]),(wx[ix+1],wy[iy+1])))
            ax.add_collection(LineCollection(edges,colors=['#252525'],linewidths=.65))
            radius = target['radius']
            ax.add_patch(Ellipse(target['pose'], 2*radius,
                2*radius*(2/3 if target['shape']=='ellipse' else 1),
                angle=target.get('orientation', 0)*180/np.pi,
                fill=False, edgecolor='#151515', lw=.8))
            ax.set_xlim(-1.5, 1.5)
            ax.set_ylim(.4, 2.8)
            ax.set_aspect('equal')
            ax.set_title(f'Texture {"AB"[tex]}: {"Re" if component==0 else "Im"} u_sc')
            ax.set_xlabel(r'$x/\lambda_0$')
            ax.set_ylabel(r'$y/\lambda_0$')
            fig.colorbar(collection, ax=ax, shrink=.8, label='Scattered field')
    table = raw_csv(name, samples)
    (DELEGATED/'P3_FIELD_RECONSTRUCTION.json').write_text(json.dumps(schemas, indent=2)+'\n')
    save(fig, name, 'Saved complex P3 scattered field: parent 0, eta=1/16, f=f0, TX0',
        sources, [table], '从保存的自由P3系数、网格及边界自由度重建局部三次多项式，在已知三角形内采样显示散射场实虚部；未将自由DOF误当网格顶点。',
        '仅可视化已保存场；非重新求解。每个细分三角形为真实P3值的平色显示，不声称连续场已认证。')


def symmetry():
    rel = 'results/symmetry/RESULTS.json'
    rows = read_json(rel)
    if not rows:
        raise FileNotFoundError(rel+' has no saved cases')
    name = 'F05_SYMMETRY_INFORMATION_AND_BREAKS'
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.7), layout='constrained')
    x = np.arange(len(rows))
    labels = [r['case'].replace('_', '\n') for r in rows]
    for offset, key, label, color in [(-.17, 'odd_wall_relative', 'Odd wall leakage', COLORS[0]),
        (.17, 'even_Bx_relative', 'Even horizontal-task leakage', COLORS[1])]:
        axes[0].bar(x+offset, [r[key] for r in rows], .34, label=label, color=color)
    axes[0].set_xticks(x, labels)
    axes[0].legend(fontsize=8)
    axes[1].bar(x, [r['target_retention'] for r in rows], color=COLORS[:len(rows)])
    axes[1].set_xticks(x, labels)
    axes[1].set_ylim(0, 1.08)
    raw = []
    for k, row in enumerate(rows):
        axes[2].plot([v['h'] for v in row['FD']], [v['relative'] for v in row['FD']],
            'o-', color=COLORS[k], label=row['case'].replace('_', ' '))
        for fd in row['FD']:
            raw.append({k:v for k,v in row.items() if k not in ('FD', 'singular_values')}
                | {'FD_h':fd['h'], 'FD_relative_error':fd['relative'],
                   'finite_odd_signal':fd['finite_odd_signal']})
    axes[0].set_title('Measured parity leakage')
    axes[1].set_title('Horizontal derivative after nuisance profiling')
    axes[2].set_title('Native derivative vs finite differences')
    decorate(axes[0], ylabel='Relative norm', log=True)
    decorate(axes[1], ylabel='Squared-norm retention')
    decorate(axes[2], 'Finite-difference half-step / lambda0', 'Relative derivative difference', log=True)
    axes[2].set_xscale('log')
    axes[2].legend(fontsize=7.5)
    table = raw_csv(name, raw)
    save(fig, name, 'Centered horizontal symmetry: available saved controls',
        [rel], [table], '显示实际奇偶泄漏、nuisance投影后横向导数保留率及导数差分审核，非对称与左右gain扰动按已保存case全部列出。',
        '中心横向局部任务控制；机器最小二乘rank不是严格rank或全目标定位证书。当前已保存case数='+str(len(rows))+'。')


def e3_distances():
    static_rel = 'results/e3_logic/offline_scores/OFFLINE_STATIC_DISTANCES.npz'
    rel = 'results/e3_logic/offline_scores/OFFLINE_DISTANCES_parent_00_repeat_00.npz'
    with np.load(ROOT/static_rel) as z:
        static = {k:z[k].copy() for k in z.files}
    with np.load(ROOT/rel) as z:
        distance = {k:z[k].copy() for k in z.files}
    name = 'F09_OLD_BANK_FOUR_DISTANCE_NOTIONS'
    metrics = [('complete_response_whitened_real_RMS', 'Full response distance', static),
        ('target_increment_whitened_real_RMS', 'Target-increment distance', static),
        ('contrast_residual_aware_oscillation', 'Residual-aware contrast oscillation', distance)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
    selected = distance['sampled_pair_mask']
    for ax, (key, label, values) in zip(axes, metrics):
        ax.scatter(values[key][~selected], distance['conditional_target_TV'][~selected],
            s=5, color=COLORS[0], alpha=.25, rasterized=True, label='All other pairs: offline')
        ax.scatter(values[key][selected], distance['conditional_target_TV'][selected],
            s=7, color=COLORS[1], alpha=.55, rasterized=True, label='Queried-wall pairs')
        ax.axhline(.02, color='#555', ls='--', lw=1)
        ax.set_xscale('symlog', linthresh=1e-3)
        ax.set_xlim(left=0.)
        ax.set_ylim(-.02, 1.04)
        decorate(ax, label, 'Conditional target posterior TV')
        ax.legend(fontsize=7, loc='best')
    raw=[]
    for i in range(len(selected)):
        raw.append(dict(parent=0, repeat=0, left_wall=int(static['pair_left'][i]),
            right_wall=int(static['pair_right'][i]), sampled_pair=bool(selected[i]),
            full_response_RMS=static['complete_response_whitened_real_RMS'][i],
            target_increment_RMS=static['target_increment_whitened_real_RMS'][i],
            background_residual_projection_norm=static['background_difference_projection_into_full_increment_space_norm'][i],
            residual_aware_contrast_oscillation=distance['contrast_residual_aware_oscillation'][i],
            conditional_target_TV=distance['conditional_target_TV'][i]))
    table=raw_csv(name, raw)
    save(fig, name, 'Old-bank diagnostic: distinct notions of background similarity',
        [static_rel,rel], [table], '固定展示parent00/repeat00的全部6903个墙pair，比较完整响应、目标增量、含残差对比与条件目标TV。',
        '完整表诊断仅在输出封存后计算，不能归为在线节省；横轴单位不同，不比较其数值大小。')


def e3_quality_cost():
    directory='results/e3_logic/offline_scores/'
    parents=read_csv(directory+'PARENT_RESULTS.csv')
    family=read_csv(directory+'FAMILY_RESULTS.csv')
    diag=read_csv(directory+'CONDITIONAL_DIAGNOSTICS.csv')
    summary=read_json(directory+'SUMMARY.json')
    methods=list(dict.fromkeys(r['method'] for r in parents))
    labels=['Response interpolation', 'Anchor-weighted contrast', 'Sampled consensus']
    name='F10_OLD_BANK_QUALITY_QUERIES_AND_FAILURES'
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for m,label,color in zip(methods,labels,COLORS):
        rr=sorted([r for r in parents if r['method']==m],key=lambda r:int(r['parent']))
        axes[0,0].plot([int(r['parent']) for r in rr],
            [float(r['target_TV_mean']) for r in rr], 'o-', color=color, ms=4,label=label)
    axes[0,0].axhline(.02,color='#555',ls='--',label='Target TV gate .02')
    axes[0,0].legend(fontsize=7.5)
    means=[next(r for r in family if r['method']==m and r['grouping']=='all') for m in methods]
    axes[0,1].bar(np.arange(3),[float(r['position_GT_RMSE']) for r in means],color=COLORS[:3])
    axes[0,1].axhline(float(means[0]['reference_GT_RMSE']),color='#555',ls='--',label='Full finite reference')
    axes[0,1].set_xticks(np.arange(3),['Response\ninterpolation','Anchor\ncontrast','Sampled\nconsensus'])
    axes[0,1].legend(fontsize=8)
    queries=[summary['shared_online_queries_per_case']]*3+[2950]
    axes[1,0].bar(np.arange(4),queries,color=COLORS)
    axes[1,0].set_xticks(np.arange(4),['Response','Contrast','Consensus','Full reference'])
    for i,v in enumerate(queries):
        axes[1,0].text(i,v+40,str(v),ha='center',fontsize=9)
    axes[1,0].text(.04,.83,'New FE / RHS = 0 / 0\nHistorical bank: 2950 FE + 11800 RHS',
        transform=axes[1,0].transAxes,fontsize=9)
    x=np.arange(len(diag))
    axes[1,1].plot(x,[float(r['all_conditional_TV_diameter']) for r in diag],
        'o',color=COLORS[0],label='Full-table conditional TV diameter')
    axes[1,1].plot(x,[float(r['unqueried_background_posterior_mass']) for r in diag],
        'x',color=COLORS[1],label='Unqueried posterior background mass')
    axes[1,1].axhline(.02,color='#555',ls='--')
    axes[1,1].set_ylim(-.04,1.07)
    axes[1,1].legend(fontsize=7.5)
    for ax,title in zip(axes.ravel(),['Parent means after two noise repeats',
        'Absolute error remains separate from regret','Logical queries per case; replay only',
        'All 24 cases: retained negative evidence']):ax.set_title(title)
    decorate(axes[0,0],'Parent index','Target TV (linear below 1e-6)')
    axes[0,0].set_yscale('symlog',linthresh=1e-6)
    axes[0,0].set_ylim(0.,1.2)
    decorate(axes[0,1],ylabel='Position RMSE / lambda0')
    decorate(axes[1,0],ylabel='Logical atom queries')
    decorate(axes[1,1],'Case index (parent-major, repeat-minor)','Probability / TV')
    tables=[raw_csv(name+'_PARENTS',parents),raw_csv(name+'_FAMILY',family),
        raw_csv(name+'_CONDITIONAL',diag),raw_csv(name+'_QUERIES',[
            dict(method=m,queries=q,new_FE=0,new_RHS=0) for m,q in zip(methods+['FULL_REFERENCE'],queries)])]
    save(fig,name,'Old-bank replay: target quality, actual queries, and consensus failure evidence',
        [directory+p for p in ('PARENT_RESULTS.csv','FAMILY_RESULTS.csv','CONDITIONAL_DIAGNOSTICS.csv','SUMMARY.json')],tables,
        '保留全部parent与两noise重复，点误差、posterior误差、查询数及未查质量分开显示；没有隐藏空墙或极差parent。',
        '查询图不是PDE速度图；真实新增FE/RHS为0。没有在线一致性证书，冷PDE支路未触发。')


def imaging():
    summary=read_json('results/imaging/RESULTS.json')
    scores=read_csv('results/imaging/TARGET_RESULTS.csv')
    with np.load(ROOT/'results/imaging/ALL_POSTERIORS_OFFLINE.npz') as z:
        recorded=json.loads(z['rows'].item())
        posterior=z['posterior'].copy()
    saved={(int(r['case']),r['contract'],r['method']):p for r,p in zip(recorded,posterior)}
    observations=np.load(ROOT/'results/imaging/OBSERVATIONS.npy')
    truth=read_json('results/imaging/TRUTH_OFFLINE.json')
    fs=(.85,1.,1.15)
    bank=np.array([[[np.load(ROOT/f'results/imaging/macro_g{g}_f{f:g}_q{j}.npz')['data']
        for f in fs] for j in range(12)] for g in range(6)])
    wall_fields=[]
    for g in range(6):
        with np.load(ROOT/f'results/imaging/macro_g{g}_f1_q0.npz') as z:
            wall_fields.append(z['eps'].copy()+1j*z['ell'].copy())
            macro_x,macro_y=z['x'].copy(),z['y'].copy()
    wall_fields=np.array(wall_fields)
    prior_mean=wall_fields.mean(axis=0)
    prior_var=float(np.mean(abs(wall_fields-prior_mean)**2))
    derived=[]
    for case,(y,tt) in enumerate(zip(observations,truth)):
        wall_rel=f'results/micro/WALL_p{tt["wall_parent"]}_n16_t{tt["texture"]}.npz'
        with np.load(ROOT/wall_rel) as z:
            true_field=z['eps'].copy()+1j*z['ell'].copy()
            xx,yy=(z['x'][:-1]+z['x'][1:])/2,(z['y'][:-1]+z['y'][1:])/2
        ix=np.searchsorted(macro_x,xx,side='right')-1
        iy=np.searchsorted(macro_y,yy,side='right')-1
        for acq in ('A0','A1','A2'):
            def select(a):
                return a[...,1,:2,:] if acq=='A0' else a[...,1,:,:] if acq=='A1' else a
            residual=select(bank)-select(y)
            ll=-.5*np.sum(abs(residual)**2,axis=tuple(range(2,residual.ndim)))/.01**2
            q=np.exp(ll-ll.max());q/=q.sum();weights=q.sum(axis=1)
            mean=np.einsum('g,gij->ij',weights,wall_fields)
            variance=float(np.mean(np.einsum('g,gij->ij',weights,abs(wall_fields-mean)**2)))
            error=float(np.sqrt(np.mean(abs(mean[iy[:,None],ix[None,:]]-true_field)**2)))
            sr=next(r for r in scores if int(r['case'])==case and r['contract']==acq
                and r['method']=='FINITE_MACRO_TARGET_MARGINAL')
            mismatch=float(np.max(abs(q.sum(axis=0)-saved[(case,acq,'FINITE_MACRO_TARGET_MARGINAL')])))
            derived.append(dict(case=case,parent=tt['parent'],noise=tt['noise'],
                wall_parent=tt['wall_parent'],texture=tt['texture'],contract=acq,
                position_MAP_error=float(sr['position_error']),
                posterior_mean_position_error=float(sr['mean_position_error']),
                wall_posterior_variance=variance,wall_prior_variance=prior_var,
                wall_variance_fraction=variance/prior_var,
                posterior_mean_wall_RMS_vs_true_micro=error,
                target_marginal_reproduction_maxabs=mismatch,
                **{f'wall_weight_{g}':float(w) for g,w in enumerate(weights)}))
    parent_rows=[]
    for acq in ('A0','A1','A2'):
        for parent in range(12):
            rr=[r for r in derived if r['contract']==acq and r['parent']==parent]
            parent_rows.append(dict(contract=acq,parent=parent,noise_repeats=len(rr),
                position_MAP_RMSE=np.sqrt(np.mean([r['position_MAP_error']**2 for r in rr])),
                posterior_mean_position_RMSE=np.sqrt(np.mean([r['posterior_mean_position_error']**2 for r in rr])),
                wall_variance_fraction_mean=np.mean([r['wall_variance_fraction'] for r in rr]),
                wall_RMS_mean=np.mean([r['posterior_mean_wall_RMS_vs_true_micro'] for r in rr])))
    name='F06_NESTED_TARGET_ERROR_AND_WALL_UNCERTAINTY'
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    methods=[('FINITE_MACRO_TARGET_MARGINAL','Finite macro marginal'),
        ('MEASUREMENT_SELECTED_MACRO','Selected macro wall'),('WALL_FIRST_FIXED','Wall first, fixed'),
        ('KNOWN_WALL_ORACLE','Known-wall oracle'),('PRIOR_ONLY','Prior only')]
    for (method,label),color in zip(methods,COLORS+['#747474']):
        rr=[next(r for r in summary['summary'] if r['contract']==acq and r['method']==method)
            for acq in ('A0','A1','A2')]
        axes[0].plot(range(3),[r['position_RMSE'] for r in rr],'o-',color=color,label=label)
    axes[0].axhline(.1,color='#555',ls='--',label='Position gate 0.1 lambda0')
    axes[0].set_xticks(range(3),['A0: 64 real','A1: 128 real','A2: 384 real'])
    axes[0].legend(fontsize=7.5)
    axes[0].set_title('Absolute MAP position error, all 24 cases')
    for acq,color,marker in zip(('A0','A1','A2'),COLORS,('o','s','^')):
        rr=[r for r in parent_rows if r['contract']==acq]
        axes[1].scatter([r['wall_variance_fraction_mean'] for r in rr],
            [r['posterior_mean_position_RMSE'] for r in rr],color=color,marker=marker,
            s=38,alpha=.8,label=acq)
    axes[1].axhline(.1,color='#555',ls='--')
    axes[1].legend(fontsize=8)
    axes[1].set_title('12 target parents; two noise repeats averaged first')
    decorate(axes[0],ylabel='Position RMSE / lambda0')
    decorate(axes[1],'Finite macro-wall variance / prior variance',
        'Posterior-mean position RMSE / lambda0')
    tables=[raw_csv(name+'_SAVED_SUMMARY',summary['summary']),raw_csv(name+'_CASES',derived),
        raw_csv(name+'_PARENTS',parent_rows)]
    save(fig,name,'Nested measurements: target error and remaining wall uncertainty',
        ['results/imaging/RESULTS.json','results/imaging/TARGET_RESULTS.csv',
         'results/imaging/ALL_POSTERIORS_OFFLINE.npz','results/imaging/OBSERVATIONS.npy',
         'results/imaging/TRUTH_OFFLINE.json','results/imaging/macro_g*_f*_q*.npz'],tables,
        '左图直接使用保存的绝对MAP误差；右图在封存后由相同有限likelihood重建墙边缘权重，先在每个target parent内平均两次noise，再比较墙方差与目标后验均值误差。',
        '12个target parent但只有2个微结构wall parent。墙方差属于6个宏观墙的近似模型，不是真实微结构posterior，也不保证模型正确或校准。墙均值RMS误差与目标边缘重现差均保存于CSV。')


def compensation():
    curves=read_json('results/compensation/CURVES.json')
    information=read_json('results/compensation/INFORMATION.json')
    validation=read_json('results/compensation/VALIDATION.json')['rows']
    name='F07_FINITE_COMPENSATION_AND_INFORMATION'
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.5),layout='constrained')
    for branch,style in ((0,'-'),(1,'--')):
        rr=sorted([r for r in curves if r['branch']==branch],key=lambda r:r['dx'])
        for acq,color in zip(('A0','A1','A2'),COLORS):
            axes[0,0].plot([r['dx'] for r in rr],[r['distance_'+acq] for r in rr],
                marker='o',ls=style,color=color,label=f'{acq}, branch {branch}')
        axes[0,1].plot([r['dx'] for r in rr],[r['material_RMS'] for r in rr],
            marker='o',ls=style,color=COLORS[branch],label=f'Branch {branch}')
    for r in validation:
        if r['f']==1.:
            axes[0,0].scatter([r['dx'],r['dx']],[r['pair_A0'],r['pair_all4']],
                c=[COLORS[0],COLORS[1]],edgecolors='black',marker='D',s=60,
                zorder=4,label='Independent refinement, f0')
    axes[0,0].axhline(1.,color='#555',ls=':',label='Pair-distance scale 1')
    axes[0,0].legend(fontsize=7,ncol=2)
    axes[0,1].legend(fontsize=8)
    task_names=['x','y','contrast','loss','radius']
    infrows=[]
    for j,task in enumerate(task_names):
        std=[np.sqrt(information['info'][acq]['prior_target_covariance'][j][j])
            for acq in ('A0','A1','A2')]
        axes[1,0].plot(range(3),std,'o-',label=task)
    for acq in ('A0','A1','A2'):
        row=information['info'][acq]
        for j,task in enumerate(task_names):
            infrows.append(dict(contract=acq,task=task,
                standardized_prior_posterior_std=np.sqrt(row['prior_target_covariance'][j][j]),
                scaled_data_target_Keffective_diagonal=row['data_target_Keffective'][j][j],
                wall_variance_retained=row['prior_wall_variance_retained'],
                float_rank=row['algebraic_float_rank'],rows=row['rows'],columns=row['columns'],
                exact_rank_status=row['exact_rank']))
    wallfraction=[information['info'][acq]['prior_wall_variance_retained'] for acq in ('A0','A1','A2')]
    axes[1,1].bar(range(3),wallfraction,color=COLORS[:3])
    for i,acq in enumerate(('A0','A1','A2')):
        rr=information['info'][acq]
        axes[1,1].text(i,wallfraction[i]+.025,f'{100*wallfraction[i]:.1f}%',ha='center')
        axes[1,1].text(i,.08,f'Float rank\n{rr["algebraic_float_rank"]} / {rr["columns"]}',
            ha='center',fontsize=8,color='white')
    for ax in axes[1]:ax.set_xticks(range(3),['A0','A1','A2'])
    axes[1,0].set_ylim(0,1.06);axes[1,1].set_ylim(0,1.08)
    axes[1,0].legend(fontsize=8,ncol=2)
    titles=['One world, two legal finite compensation branches',
        'Actual material changes of the saved endpoints',
        'Proper tangent prior: target standard deviation fraction',
        'Proper tangent prior: wall variance still retained']
    for ax,title in zip(axes.ravel(),titles):ax.set_title(title)
    decorate(axes[0,0],'Target horizontal displacement / lambda0','Whitened pair distance',log=True)
    decorate(axes[0,1],'Target horizontal displacement / lambda0','Wall material RMS change')
    decorate(axes[1,0],ylabel='Standard deviation / tangent prior std')
    decorate(axes[1,1],ylabel='Wall variance / tangent prior variance')
    tables=[raw_csv(name+'_CURVES',[{k:v for k,v in r.items() if k!='history'} for r in curves]),
        raw_csv(name+'_INFORMATION',infrows),raw_csv(name+'_VALIDATION',validation),
        raw_csv(name+'_DERIVATIVE_AUDIT',information['FD'])]
    save(fig,name,'Finite full-wave compensation and genuinely nested observations',
        ['results/compensation/CURVES.json','results/compensation/INFORMATION.json',
         'results/compensation/VALIDATION.json'],tables,
        '一个基准world的两个合法延拓分支完整显示；新增观测在同一端点上评分。独立网格f0复核点另用黑边菱形显示。下排为独立proper切线prior诊断。',
        '两个分支不是两个独立物理parent。找到的端点仅给歧义存在证据，不给全局profile最优或上界；float rank不是认证精确rank，切线Gaussian不是boxed nonlinear posterior。')


def return_operator():
    original=read_csv('results/e3_return_operator/PER_RANK_ERRORS.csv')
    repaired=read_csv('results/e3_return_transpose_repair/PER_RANK_ERRORS.csv')
    name='F08_RETURN_RANK_ERROR_SPECTRUM_AND_PAID_COST'
    fig,axes=plt.subplots(2,2,figsize=(13,8.5),layout='constrained')
    series=[(original,'return_only_exact_target_self','Original adjoint projection',COLORS[1],'--'),
        (repaired,'return_transpose_exact_target_self','Transpose repair, exact target self',COLORS[0],'-'),
        (repaired,'whole_C_transpose_same_return_basis_diagnostic','Whole C, same H-fitted basis',COLORS[2],':')]
    aggregates=[]
    for table,method,label,color,style in series:
        values=[]
        for rank in (4,8,16,32):
            rr=[r for r in table if r['method']==method and int(r['rank'])==rank]
            row=dict(method=method,rank=rank,cases=len(rr),
                max_white=max(float(r['white_response_norm']) for r in rr),
                max_induced_ROI_L2=max(float(r['ROI_target_induced_L2_relative']) for r in rr),
                max_true_current_closure=max(float(r['true_feedback_current_relative']) for r in rr))
            aggregates.append(row);values.append(row)
        axes[0,0].plot([r['rank'] for r in values],[r['max_white'] for r in values],
            marker='o',ls=style,color=color,label=label)
        if method!='whole_C_transpose_same_return_basis_diagnostic':
            axes[0,1].plot([r['rank'] for r in values],[r['max_induced_ROI_L2'] for r in values],
                marker='o',ls=style,color=color,label=label+': ROI field')
            axes[0,1].plot([r['rank'] for r in values],[r['max_true_current_closure'] for r in values],
                marker='x',ls=style,color=color,alpha=.7,label=label+': closure')
    axes[0,0].axhline(.1,color='#555',ls='--',label='Response qualification scale .1')
    for ax in axes[0]:
        ax.set_xticks([4,8,16,32]);ax.legend(fontsize=7)
    spectra=[]
    for wall,style in ((0,'-'),(2,'--')):
        rel=f'results/e3_return_transpose_repair/wall_{wall}/RAW_TRANSPOSE_ACTIONS.npz'
        with np.load(ROOT/rel) as z:
            for key,label,color in [('singular_values_Hcore','H projected core',COLORS[0]),
                ('singular_values_Ccore','Cw projected core',COLORS[1])]:
                s=z[key].copy()
                axes[1,0].plot(np.arange(1,len(s)+1),s/s[0],ls=style,color=color,
                    label=f'Wall {wall}: {label}')
                spectra += [dict(wall=wall,core=label,mode=i+1,singular_value=float(v),
                    normalized_by_own_first=float(v/s[0]),source=rel) for i,v in enumerate(s)]
    cost=[]
    for label,rel in [('Original task','results/jobs/return-operator-frozen/RECEIPT.json'),
        ('Transpose repair','results/jobs/return-transpose-repair-frozen/RECEIPT.json')]:
        r=read_json(rel)
        rhs,fe=(708,19) if label=='Original task' else (256,3)
        cost.append(dict(component=label,RHS=rhs,FE=fe,wrapper_CPU=r['process_cpu_seconds'],source=rel))
    workload=[dict(component='Original task',RHS=708,FE=19),
        dict(component='Transpose repair',RHS=256,FE=3),
        dict(component='Task total',RHS=964,FE=22),
        dict(component='8 direct references',RHS=16,FE=8)]
    axes[1,1].bar(range(4),[r['RHS'] for r in workload],color=COLORS)
    axes[1,1].set_xticks(range(4),['Original\ntask','Transpose\nrepair','Task\ntotal','8 direct\nreferences'])
    for i,r in enumerate(workload):
        axes[1,1].text(i,r['RHS']+20,f'{r["RHS"]} RHS\n{r["FE"]} FE',ha='center',fontsize=8)
    axes[1,1].set_ylim(0,1110)
    for ax,title in zip(axes.ravel(),['Worst of all 8 wall-target cases',
        'Full-feedback field and current defects',
        '32 x 32 projected cores: coefficient metric',
        'Paid workloads; separate scientific purposes']):ax.set_title(title)
    decorate(axes[0,0],'Frozen nested rank','Max whitened response norm',log=True)
    decorate(axes[0,1],'Frozen nested rank','Max relative error / defect',log=True)
    decorate(axes[1,0],'Projected-core singular mode','Singular value / own first value',log=True)
    decorate(axes[1,1],ylabel='Complex RHS requests')
    axes[1,0].legend(fontsize=7)
    tables=[raw_csv(name+'_ORIGINAL',original),raw_csv(name+'_REPAIRED',repaired),
        raw_csv(name+'_AGGREGATES',aggregates),raw_csv(name+'_CORE_SPECTRA',spectra),
        raw_csv(name+'_WORKLOAD',workload),raw_csv(name+'_WRAPPER_COST',cost)]
    save(fig,name,'Return representation: preserved first attempt and one transpose repair',
        ['results/e3_return_operator/PER_RANK_ERRORS.csv','results/e3_return_transpose_repair/PER_RANK_ERRORS.csv',
         'results/e3_return_transpose_repair/wall_*/RAW_TRANSPOSE_ACTIONS.npz',
         'results/jobs/return-operator-frozen/RECEIPT.json','results/jobs/return-transpose-repair-frozen/RECEIPT.json'],tables,
        '2墙×4目标，共8组合；2train和2heldout全部保留。原始负面平台和一次transpose修复同图显示，各秩不挑选最佳case。谱仅来自保存的小投影core。',
        'whole-C使用同一个仅为H训练的基底，不是最优whole-C对照。Euclidean FE系数core谱不是完整H或物理L2谱；没有operator norm、posterior、连续PDE证书或等成本加速结论。原始/修复/总任务含setup与validation，direct reference只承担评分。')


def feedback():
    name='F12_FULL_SELF_AND_BORN_FROM_PAID_CACHE'
    rows=[];receivers=[];sources=[]
    for wall in (0,2):
        base=ROOT/f'results/e3_return_operator/wall_{wall}'
        with np.load(base/'BACKGROUND.npz') as z:
            background=z['d_wall'].copy();Bw=z['Bw'].copy()
        sources.append(str((base/'BACKGROUND.npz').relative_to(ROOT)))
        for directory in sorted(base.glob('target_*')):
            target=int(directory.name.split('_')[1])
            with np.load(directory/'EXACT_SELF_TEMPLATE.npz') as z:
                q,j0=z['q'].copy(),z['j0'].copy()
            with np.load(directory/'DIRECT_REFERENCE.npz') as z:
                full=z['data'].copy();current=z['j'].copy()
            born=background+(Bw@q).T
            self_only=background+(Bw@j0).T
            denominator=np.linalg.norm(full-background)
            if not denominator>0:raise ValueError('vanishing full target-induced response')
            rows.append(dict(wall=wall,target=target,
                born_induced_relative_error=np.linalg.norm(born-full)/denominator,
                self_only_induced_relative_error=np.linalg.norm(self_only-full)/denominator,
                born_white_error=np.linalg.norm(born-full)/.01,
                self_only_white_error=np.linalg.norm(self_only-full)/.01,
                full_induced_response_norm=denominator,
                receiver_reconstruction_maxabs=np.max(abs(background+(Bw@current).T-full))))
            if wall==0 and target==1:
                for rx in range(full.shape[1]):
                    receivers.append(dict(wall=wall,target=target,source=0,receiver=rx,
                        full_induced_real=(full-background)[0,rx].real,
                        full_induced_imag=(full-background)[0,rx].imag,
                        born_induced_real=(born-background)[0,rx].real,
                        born_induced_imag=(born-background)[0,rx].imag,
                        exact_self_only_real=(self_only-background)[0,rx].real,
                        exact_self_only_imag=(self_only-background)[0,rx].imag))
    fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    x=np.arange(len(rows))
    axes[0].bar(x-.17,[r['born_induced_relative_error'] for r in rows],.34,
        color=COLORS[1],label='Wall-exact Born insertion')
    axes[0].bar(x+.17,[r['self_only_induced_relative_error'] for r in rows],.34,
        color=COLORS[0],label='Exact target self; no wall return')
    axes[0].set_xticks(x,[f'W{r["wall"]}\nT{r["target"]}' for r in rows])
    axes[0].legend(fontsize=7.5)
    for column,component in ((1,'real'),(2,'imag')):
        for prefix,label,color,style in [('full_induced','Full feedback',COLORS[0],'-'),
            ('born_induced','Born insertion',COLORS[1],'--'),
            ('exact_self_only','Exact target self only',COLORS[2],':')]:
            axes[column].plot([r['receiver'] for r in receivers],
                [r[prefix+'_'+component] for r in receivers],ls=style,color=color,label=label)
        axes[column].legend(fontsize=7.5)
        axes[column].set_title(f'W0 / T1 strong ellipse, TX0: {component}')
        decorate(axes[column],'Receiver index',f'{component.capitalize()} target-induced response')
    axes[0].set_title('All 8 paid wall-target configurations')
    decorate(axes[0],ylabel='Error / full target-induced response norm',log=True)
    tables=[raw_csv(name+'_ALL_CASES',rows),raw_csv(name+'_RECEIVERS',receivers)]
    save(fig,name,'Born, exact target self scattering, and complete wall-return feedback',
        sources+['results/e3_return_operator/wall_*/target_*/EXACT_SELF_TEMPLATE.npz',
                 'results/e3_return_operator/wall_*/target_*/DIRECT_REFERENCE.npz'],tables,
        '只重组已付费Bw、q、j0与full参考：Born为dW+Bw q，自散射控制为dW+Bw j0，full来自保存的完整FEM。固定展示W0/T1强椭圆的真实接收响应实虚部。',
        '零新增PDE。2墙×4目标的有限共同空间控制；用完整目标诱导响应作分母，非绝对总场；不是均匀物理或posterior证书。')


def failure_boundary():
    curves=read_json('results/compensation/CURVES.json')
    symmetry_rows=read_json('results/symmetry/RESULTS.json')
    maximum=max(r['dx'] for r in curves)
    endpoints=[r for r in curves if r['dx']==maximum]
    name='F11_MEASURED_FAILURE_BOUNDARIES'
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
    x=np.arange(3)
    for row,offset,color in zip(endpoints,(-.17,.17),COLORS):
        values=[row['distance_'+a] for a in ('A0','A1','A2')]
        axes[0].bar(x+offset,values,.34,color=color,label=f'Compensation branch {row["branch"]}')
    axes[0].axhline(1.,color='#555',ls='--',label='Pair-distance scale 1')
    axes[0].set_xticks(x,['A0','A1','A2']);axes[0].legend(fontsize=8)
    axes[0].set_title(f'Same legal endpoint at dx={maximum:g} lambda0')
    axes[1].bar(range(len(symmetry_rows)),[1-r['target_retention'] for r in symmetry_rows],
        color=COLORS[:len(symmetry_rows)])
    axes[1].set_xticks(range(len(symmetry_rows)),[r['case'].replace('_','\n') for r in symmetry_rows])
    axes[1].set_title('Fraction of horizontal derivative energy lost')
    decorate(axes[0],ylabel='Whitened pair distance',log=True)
    decorate(axes[1],ylabel='1 - target retention')
    raw=[dict(kind='compensation',branch=r['branch'],dx=r['dx'],
        A0=r['distance_A0'],A1=r['distance_A1'],A2=r['distance_A2']) for r in endpoints]
    raw += [dict(kind='parity',case=r['case'],retention=r['target_retention'],
        lost_fraction=1-r['target_retention'],odd_wall_relative=r['odd_wall_relative'],
        even_Bx_relative=r['even_Bx_relative']) for r in symmetry_rows]
    table=raw_csv(name,raw)
    save(fig,name,'Observed boundaries of narrow measurement equivalence and parity protection',
        ['results/compensation/CURVES.json','results/symmetry/RESULTS.json'],[table],
        '同一个实际补偿端点在新增孔径/频率下不再保持原A0噪声尺度等效；非对称墙或差分gain可削弱中心横向任务保护。',
        '这里按同一已找到端点比较，不能断言A1/A2下不存在另一补偿端点。对称压力仅一个物理world的三个控制，不是三个独立重复或一般unknown-wall失败定理。')


def micro_finite():
    summary=read_json('results/e2_micro_reference/SUMMARY.json')
    cases=read_csv('results/e2_micro_reference/CASE_RESULTS.csv')
    consensus=read_json('results/e2_micro_reference/MICRO_CONSENSUS_OFFLINE.json')
    rows=[next(r for r in summary['summary'] if r['contract']==a) for a in ('A0','A1','A2')]
    tv=[next(r for r in consensus['summary'] if r['contract']==a) for a in ('A0','A1','A2')]
    name='F13_CONTROLLED_MICRO_RISK_TEXTURE_AND_TV_TAILS'
    fig,axes=plt.subplots(2,3,figsize=(15,8.6),layout='constrained')
    x=np.arange(3)
    for key,label,color in [('position_MAP_RMSE','MAP position RMSE',COLORS[0]),
        ('position_mean_RMSE','Posterior-mean position RMSE',COLORS[1])]:
        axes[0,0].plot(x,[r[key] for r in rows],'o-',color=color,label=label)
    axes[0,0].axhline(.1,color='#555',ls='--',label='Position gate .1 lambda0')
    axes[0,1].plot(x,[r['posterior_position_risk'] for r in rows],'o-',color=COLORS[0],label='Posterior position risk')
    axes[0,1].plot(x,[r['prior_position_risk'] for r in rows],'--',color='#777',label='Prior position risk')
    axes[0,2].plot(x,[r['texture_entropy_bits'] for r in rows],'o-',color=COLORS[0],label='Texture entropy')
    axes[0,2].plot(x,[r['conditional_texture_entropy_true_macro_bits'] for r in rows],
        's--',color=COLORS[1],label='Given true macro: OFFLINE')
    axes[0,2].set_ylim(0,1.05)
    axes[1,0].plot(x,[r['global_wall_pixel_variance_ratio_to_prior'] for r in rows],
        'o-',color=COLORS[0],label='Global 4-wall finite posterior')
    axes[1,0].plot(x,[r['truth_macro_pixel_variance_ratio_to_prior'] for r in rows],
        's--',color=COLORS[1],label='Given true macro: OFFLINE')
    axes[1,0].set_ylim(0,1.06)
    tvrows=[]
    for key,label,color,style in [('within_macro_conditional_TV','All within-macro pairs (48)',COLORS[0],'-'),
        ('true_macro_conditional_TV','True-macro pair only: OFFLINE (24)',COLORS[1],'--'),
        ('all_four_background_conditional_diameter','All four backgrounds (24)',COLORS[2],':')]:
        axes[1,1].plot(x,[r[key]['maximum'] for r in tv],marker='o',ls=style,color=color,label=label)
    axes[1,1].axhline(.02,color='#555',ls='--',label='Target TV scale .02')
    for r in tv:
        tvrows.append(dict(contract=r['contract'],within_macro_pair_observations=r['within_macro_pair_observations'],
            all_within_macro_TV_max=r['within_macro_conditional_TV']['maximum'],
            true_macro_TV_max_OFFLINE=r['true_macro_conditional_TV']['maximum'],
            all_four_background_TV_diameter_max=r['all_four_background_conditional_diameter']['maximum'],
            all_within_macro_TV_above_002_count=r['within_macro_TV_above_002_count'],
            true_macro_TV_above_002_count_OFFLINE=r['true_macro_TV_above_002_count']))
    tails=consensus['all_failures_above_002']
    for i,r in enumerate(tails):
        axes[1,2].bar(i,r['conditional_target_TV'],color=COLORS[('A0','A1','A2').index(r['contract'])])
        axes[1,2].text(i,r['conditional_target_TV']+.0008,f'{r["conditional_target_TV"]:.5f}',ha='center',fontsize=7)
    axes[1,2].set_xticks(range(len(tails)),
        [f'{r["contract"]}\nP{r["parent"]} n{r["noise"]}\nmacro {r["macro_a"]}' for r in tails])
    axes[1,2].axhline(.02,color='#555',ls='--')
    axes[1,2].set_ylim(0,.044)
    for ax in list(axes[0])+list(axes[1,:2]):ax.set_xticks(x,['A0','A1','A2'])
    titles=['Absolute target error, 24 observations', 'Target posterior risk, same finite prior',
        'Texture alternatives remain unresolved', 'Actual wall pixel variance still retained',
        'Maxima over every recorded micro group', 'All five within-macro TV tails above .02']
    for ax,title in zip(axes.ravel(),titles):ax.set_title(title)
    decorate(axes[0,0],ylabel='Position RMSE / lambda0')
    decorate(axes[0,1],ylabel='Position risk / lambda0 squared')
    decorate(axes[0,2],ylabel='Texture entropy / bits')
    decorate(axes[1,0],ylabel='Pixel variance / corresponding prior variance')
    decorate(axes[1,1],ylabel='Conditional target TV maximum',log=True)
    decorate(axes[1,2],ylabel='Conditional target TV')
    for ax in list(axes[0])+list(axes[1,:2]):ax.legend(fontsize=7)
    tables=[raw_csv(name+'_SUMMARY',rows),raw_csv(name+'_CASES',cases),
        raw_csv(name+'_TV_MAXIMA',tvrows),raw_csv(name+'_ALL_POSITIVE_TAILS',tails)]
    save(fig,name,'OFFLINE controlled 4-micro-wall x 12-target finite reference: risk and texture uncertainty',
        ['results/e2_micro_reference/SUMMARY.json','results/e2_micro_reference/CASE_RESULTS.csv',
         'results/e2_micro_reference/MICRO_CONSENSUS_OFFLINE.json'],tables,
        '4个真实已付费微结构墙×12个联合目标标签的离线有限参考；12个target scenario、2个固定macro parent、每scenario两noise，共24观测。所有微结构组的TV最大值及全部5个>.02尾部逐项保留。',
        '48原子proper有限prior、同生成家族控制；没有修改在线prior。给定真实macro的熵/方差/TV都明确标OFFLINE，不是在线真值选组。全表均已查询，不能据此声称在线停止、查询节省、连续posterior或统一Maxwell证书；两种texture的损耗相同，不能由此宣称loss不确定。')


def final_compensation():
    coarse=read_json('results/compensation/CURVES.json')
    extended=read_json('results/compensation_extended/RESULTS.json')['curves']
    final=read_json('results/compensation_final/RESULTS.json')
    third=read_json('results/final_pair_refinement/RESULTS.json')
    baseline_rel='results/compensation/BASE.npz'
    with np.load(ROOT/baseline_rel) as z:
        base_data=z['data'].copy();eps0=z['eps'].copy();ell0=z['ell'].copy();wx=z['x'].copy();wy=z['y'].copy()
    with np.load(ROOT/'results/compensation_final/ENDPOINT.npz') as z:
        endpoint=z['data'].copy();eps1=z['eps'].copy();ell1=z['ell'].copy()
    def distances(base,end):
        delta=end-base
        return [float(np.linalg.norm(delta[1,:2])/.01),
            float(np.linalg.norm(delta[1])/.01),float(np.linalg.norm(delta)/.01)]
    a0=next(r for r in coarse if r['branch']==0 and r['dx']==.05)
    stages=[dict(stage='A0 optimized',**{a:a0['distance_'+a] for a in ('A0','A1','A2')},material_RMS=a0['material_RMS'])]
    for acq in ('A1','A2'):
        r=next(r for r in extended if r['contract']==acq)
        stages.append(dict(stage=acq+' reoptimized, 5 steps',**{a:r['d_'+a] for a in ('A0','A1','A2')},material_RMS=r['material_RMS']))
    stages.append(dict(stage='A2, 7 steps: primary grid',**dict(zip(('A0','A1','A2'),distances(base_data,endpoint))),material_RMS=final['material_RMS']))
    for label,folder,template in [('A2, 7 steps: refinement 2','compensation_final','INDEPENDENT_f{f}.npz'),
        ('A2, 7 steps: refinement 3','final_pair_refinement','f{f}.npz')]:
        pairs=[]
        for f in (.85,1.,1.15):
            with np.load(ROOT/f'results/{folder}'/template.format(f=f)) as z:
                pairs.append((z['base'].copy(),z['endpoint'].copy()))
        d=distances(np.array([p[0] for p in pairs]),np.array([p[1] for p in pairs]))
        stages.append(dict(stage=label,**dict(zip(('A0','A1','A2'),d)),material_RMS=final['material_RMS']))
    accepted=[];trials=[]
    five=next(r for r in extended if r['contract']=='A2')['history']
    accepted.append(dict(step=0,distance=five[0]['d'],stage='start from A0 endpoint'))
    for step,h in enumerate(five,1):
        for t in h['trials']:trials.append(dict(step=step,alpha=t['alpha'],distance=t['d']))
        accepted.append(dict(step=step,distance=min([h['d']]+[t['d'] for t in h['trials']]),stage='original five-step A2 profile'))
    for step,h in enumerate(final['history'],6):
        for t in h['trials']:trials.append(dict(step=step,alpha=t['alpha'],distance=t['d']))
        accepted.append(dict(step=step,distance=min([h['start']]+[t['d'] for t in h['trials']]),stage='two additional authorized steps'))
    name='F14_FINAL_FINITE_COMPENSATION_AND_THIRD_GRID'
    fig,axes=plt.subplots(2,3,figsize=(15,8),gridspec_kw={'height_ratios':[1.45,.65]},layout='constrained')
    for acq,color in zip(('A0','A1','A2'),COLORS):
        axes[0,0].plot(range(len(stages)),[r[acq] for r in stages],'o-',color=color,label=acq)
    axes[0,0].set_xticks(range(len(stages)),['A0 fit','A1 fit\n5 steps','A2 fit\n5 steps','A2 fit\n7 steps','Grid 2','Grid 3'])
    axes[0,0].axhline(1.,color='#555',ls='--',label='Pair scale 1')
    axes[0,0].legend(fontsize=7)
    axes[0,1].plot([r['step'] for r in accepted],[r['distance'] for r in accepted],
        'o-',color=COLORS[0],label='Accepted endpoint sequence')
    axes[0,1].scatter([r['step'] for r in trials],[r['distance'] for r in trials],
        facecolors='none',edgecolors=COLORS[1],s=35,label='Every evaluated trial')
    axes[0,1].axhline(1.,color='#555',ls='--');axes[0,1].axvline(5.5,color='#999',ls=':')
    axes[0,1].set_xticks(range(8));axes[0,1].legend(fontsize=7)
    grid_values=[stages[i]['A2'] for i in (3,4,5)]
    axes[0,2].bar(range(3),grid_values,color=COLORS[:3])
    for i,v in enumerate(grid_values):axes[0,2].text(i,v+.012,f'{v:.6f}',ha='center',fontsize=8)
    axes[0,2].set_xticks(range(3),['Primary','Refinement 2','Refinement 3'])
    axes[0,2].axhline(1.,color='#555',ls='--');axes[0,2].set_ylim(0,1.08)
    delta=eps1-eps0;bound=float(np.max(abs(delta)))
    maps=[(eps0,'Baseline wall', 'viridis',1.5,6.),
        (eps1,'Actual 7-step endpoint wall','viridis',1.5,6.),
        (delta,'Endpoint minus baseline','RdBu_r',-bound,bound)]
    material=[]
    for ax,(value,title,cmap,low,high) in zip(axes[1],maps):
        art=ax.pcolormesh(wx,wy,value,cmap=cmap,vmin=low,vmax=high,rasterized=True)
        ax.set_aspect('equal');ax.set_xlabel('x / lambda0');ax.set_ylabel('y / lambda0');ax.set_title(title)
        fig.colorbar(art,ax=ax,shrink=.8,label='epsilon_inf' if cmap=='viridis' else 'Delta epsilon_inf')
    for iy,yc in enumerate((wy[:-1]+wy[1:])/2):
        for ix,xc in enumerate((wx[:-1]+wx[1:])/2):
            material.append(dict(x=xc,y=yc,eps_base=eps0[iy,ix],eps_endpoint=eps1[iy,ix],
                ell_base=ell0[iy,ix],ell_endpoint=ell1[iy,ix]))
    axes[0,0].set_title('Each actual endpoint, scored on all acquisitions')
    axes[0,1].set_title('A2 profile: five steps, then two more')
    axes[0,2].set_title('Same final pair: all 384 real coordinates')
    decorate(axes[0,0],ylabel='Whitened pair distance',log=True)
    decorate(axes[0,1],'Bounded corrector step','A2 whitened pair distance')
    decorate(axes[0,2],ylabel='A2 whitened pair distance')
    tables=[raw_csv(name+'_ENDPOINT_STAGES',stages),raw_csv(name+'_ACCEPTED_STEPS',accepted),
        raw_csv(name+'_ALL_TRIALS',trials),raw_csv(name+'_FINAL_MATERIAL',material),
        raw_csv(name+'_THIRD_GRID',third['rows']),raw_csv(name+'_SECOND_GRID',final['rows'])]
    save(fig,name,'Actual finite compensation at dx=0.05 lambda0: reoptimization and three grids',
        ['results/compensation/CURVES.json','results/compensation/BASE.npz',
         'results/compensation_extended/RESULTS.json','results/compensation_final/RESULTS.json',
         'results/compensation_final/ENDPOINT.npz','results/compensation_final/INDEPENDENT_f*.npz',
         'results/final_pair_refinement/RESULTS.json','results/final_pair_refinement/f*.npz'],tables,
        '保留A0优化端点、A1重新优化、A2五步及七步端点；三套网格的最终同一pair直接比较。全部真实试探点、最后材料及损耗逐格导出，目标实际位移.05λ0。',
        '这是有界development profile及离散独立细化证据，不是全局profile infimum、精确连续等值族或连续误差上界。原F07/F11只说明原端点在新增测量下的变化；本图显示重新优化可重新找到噪声尺度pair，不能把旧图升级成唯一性。')


def continuous_starts():
    original_rel='results/continuous_fwi/ONLINE_OUTPUTS.json'
    third_rel='results/continuous_fwi_third/ONLINE_OUTPUTS.json'
    original=read_json(original_rel);third=read_json(third_rel)
    true_pose=np.asarray(read_json('results/continuous_fwi/TRUTH_OFFLINE.json')['target']['pose'])
    runs=[];evaluations=[]
    for rel,table,is_third in [(original_rel,original,False),(third_rel,third,True)]:
        for row in table:
            start=2 if is_third else int(row['start'])
            error=float(np.linalg.norm(np.asarray(row['best_target']['pose'])-true_pose))
            runs.append(dict(method=row['method'],public_start=start,source_start_id=row['start'],
                shape=row['shape'],best_objective=row['best_objective'],position_error=error,
                best_x=row['best_target']['pose'][0],best_y=row['best_target']['pose'][1],
                best_orientation=row['best_target']['orientation'],calls=row['calls'],stop=row['stop'],
                start_source=rel,selected_by='lowest objective within method; truth not used'))
            for h in row['history']:
                evaluations.append(dict(method=row['method'],public_start=start,call=h['call'],
                    objective=h['objective'],position_error=float(np.linalg.norm(np.asarray(h['target']['pose'])-true_pose)),
                    x=h['target']['pose'][0],y=h['target']['pose'][1],radius=h['target']['radius'],
                    contrast=h['target']['contrast'],loss=h['target']['loss'],
                    orientation=h['target']['orientation'],source=rel))
    methods=['KNOWN_WALL_ORACLE','JOINT_MACRO_CELL_FWI','EFFECTIVE_MACRO_FWI']
    chosen={method:min([r for r in runs if r['method']==method],key=lambda r:r['best_objective'])['public_start'] for method in methods}
    for r in runs:r['selected_by_objective_among_three']=r['public_start']==chosen[r['method']]
    name='F15_CONTINUOUS_ALL_STARTS_OBJECTIVE_AND_TRUE_ERROR'
    fig,axes=plt.subplots(1,3,figsize=(15,5.1),layout='constrained',sharey=True)
    labels=['Known-wall oracle','Ordinary joint macro-cell FWI','Effective macro FWI']
    for ax,method,title in zip(axes,methods,labels):
        for start,marker,color in zip((0,1,2),('o','s','^'),COLORS):
            ee=[r for r in evaluations if r['method']==method and r['public_start']==start]
            ax.plot([r['objective'] for r in ee],[r['position_error'] for r in ee],
                marker=marker,color=color,alpha=.65,lw=.8,ms=3.5,label=f'Public start S{start}: all calls')
            best=next(r for r in runs if r['method']==method and r['public_start']==start)
            ax.scatter(best['best_objective'],best['position_error'],s=65,
                marker=marker,facecolors='white',edgecolors=color,linewidth=1.7,zorder=4)
        best=next(r for r in runs if r['method']==method and r['selected_by_objective_among_three'])
        ax.scatter(best['best_objective'],best['position_error'],s=190,marker='*',
            color='#111',edgecolors='white',linewidth=.6,zorder=5,label='Selected by lowest objective')
        ax.axhline(.1,color='#555',ls='--',label='Position gate .1 lambda0')
        ax.text(.04,.97,f'Objective selects S{chosen[method]}\nAll three runs: NONCONVERGED',
            transform=ax.transAxes,va='top',fontsize=8,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
        ax.set_title(title);ax.set_ylim(0,1.62)
        decorate(ax,'Whitened least-squares objective','Actual position error / lambda0')
        ax.legend(fontsize=7,loc='lower right')
    tables=[raw_csv(name+'_ALL_RUNS',runs),raw_csv(name+'_ALL_OBJECTIVE_CALLS',evaluations)]
    save(fig,name,'Continuous development pilot: all three public starts, objective-only selection',
        [original_rel,third_rel,'results/continuous_fwi/TRUTH_OFFLINE.json',
         'results/continuous_fwi_third/PLAN.json'],tables,
        '一个离网格目标world、三个方法各三个公共起点，共9个有界run和90次真实目标函数调用。星号严格按每方法最低best objective选择；真值只用于离线误差评分。普通joint的最低objective分支位置错误，新增公共分支虽位置较好也未被该规则选中。',
        '全部NONCONVERGED_CALL_CAP，不能升级为充分收敛的成熟FWI或普遍成像失败。S2为后登记的中性ROI中心及此前缺失的允许椭圆方向，development adaptation和成本均保留；没有用GT挑起点或挑winner。')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--sections',default='micro,symmetry,e3,imaging,compensation,return_operator,feedback,failure_boundary,micro_finite,final_compensation,continuous_starts')
    parser.add_argument('--allow-missing',action='store_true')
    args=parser.parse_args()
    FIG.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    DELEGATED.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
        'axes.titlesize':10,'axes.labelsize':9,'legend.frameon':False,
        'pdf.fonttype':42,'ps.fonttype':42,'savefig.pad_inches':.12})
    sections=args.sections.split(',')
    failed=[]
    available={}
    tasks={'micro':lambda:(micro_materials(read_json('results/micro/RESULTS.json')['rows']),
        micro_scale(read_json('results/micro/RESULTS.json')['rows']),micro_validation(),micro_fields()),
        'symmetry':symmetry,'e3':lambda:(e3_distances(),e3_quality_cost()),
        'imaging':imaging,'compensation':compensation,'return_operator':return_operator,
        'feedback':feedback,'failure_boundary':failure_boundary,
        'micro_finite':micro_finite,'final_compensation':final_compensation,
        'continuous_starts':continuous_starts}
    for name in sections:
        try:
            tasks[name]()
            available[name]='PLOTTED_AVAILABLE_SAVED_RESULTS'
        except FileNotFoundError as error:
            if not args.allow_missing:
                raise
            available[name]='RESULTS_PENDING_NOT_PLOTTED'
            failed.append(dict(section=name,error=str(error)))
    manifest_path=FIG/'FIGURE_MANIFEST.json'
    previous=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    entries=previous.get('figures',{})|REGISTRY
    section_status=previous.get('sections',{})|available
    required=['micro','symmetry','e3','imaging','compensation','return_operator','feedback','failure_boundary',
        'micro_finite','final_compensation','continuous_starts']
    manifest=dict(schema='R2_SAVED_EVIDENCE_FIGURES_V1',updated_unix=time.time(),
        figures=entries,sections=previous.get('sections',{})|available,
        pending_sections=[k for k in required if section_status.get(k)!='PLOTTED_AVAILABLE_SAVED_RESULTS'],
        missing_results=failed,new_PDE=False,hash_checks=False,publication=False,
        validation='PNG/PDF saved; raw CSV exported; scientific interpretation belongs to parent review',
        rendering_note_zh='图内为英文；中文说明、数据来源和范围限制在各figure元数据。')
    manifest_path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    lines=['# R2图表清单','', '只读取保存的实际结果。没有新增PDE、摘要检查或发布。','']
    for name,item in entries.items():
        lines += [f'## {name}', '',item['explanation_zh'],'',item['scope_zh'],'',
                  'PNG: `'+item['png']+'`；PDF: `'+item['pdf']+'`。','']
    (DELEGATED/'FIGURES_ZH.md').write_text('\n'.join(lines)+'\n')
    summary=dict(figures_written_this_run=list(REGISTRY),total_figures=len(entries),
        pending_sections=manifest['pending_sections'],missing_results=failed,no_PDE=True)
    (DELEGATED/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':
    main()
