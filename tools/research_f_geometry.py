"""F0: decoder-only geometry and cached J algebra; never imports a solver.

Existing Maxwell Hessian vectors are optional. Historical scalar norms are
reported as such and never promoted to vectors or recovered with cached LU.
"""
from __future__ import annotations
import argparse, ast, csv, json, sys, time, zipfile
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.special import ndtr, ndtri

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from prior_constrained import prior  # AST-only rasterizer; no physics import

COUNTS = Counter()
STEPS = (1e-3, 5e-4, 2.5e-4)
RANK_RTOL = 1e-7


def plain(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, dict): return {k: plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [plain(v) for v in x]
    return x


def rel(x, y): return float(np.linalg.norm(x-y)/max(np.linalg.norm(y), 1e-300))
def norm(x): return float(np.linalg.norm(x))
def mass(x): return np.r_[x.real, x.imag]/24
def cosine(x, y): return float(x@y/(np.linalg.norm(x)*np.linalg.norm(y)))


def load_chart_class():
    # Load only the existing decoder class, excluding all orchestration imports.
    src = ROOT / 'src/lxz/local2d.py'
    tree = ast.parse(src.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'MaterialPrior')
    ns = dict(np=np, ndtr=ndtr, ndtri=ndtri, prior=prior, copy=__import__('copy'))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[cls], type_ignores=[])), str(src), 'exec'), ns)
    return ns['MaterialPrior']


def evaluate(chart, u):
    COUNTS['decoder_evaluations'] += 1
    x, row = chart.decode(u)
    vals = np.asarray(prior.rho_field(x, prior.background(chart.bg)))
    a = ndtr(u)
    return x, dict(active_max_cell=int(np.argmax(vals)),
        normalized_peak_gap=float(np.sort(vals)[-1]-np.sort(vals)[-2]),
        cdf_min=float(a.min()), cdf_max=float(a.max()),
        rounded_cdf_saturation=bool(np.any((a == 0) | (a == 1))))


def analytic_second(chart, u, v):
    """Exact fixed-maximum-cell directional D2 of the two-Gaussian map."""
    COUNTS['analytic_decoder_second_derivatives'] += 1
    d = chart.row['geometry']
    p = prior.geometry_points()-np.asarray(d['center'])
    co, si = np.cos(d['rotation']), np.sin(d['rotation'])
    pts = np.c_[(co*p[:, 0]+si*p[:, 1])/d['scale'],
                 (-si*p[:, 0]+co*p[:, 1])/(d['scale']*d['axis_ratio'])]
    span = chart.hi-chart.lo
    a = chart.lo+span*ndtr(u)
    ap = span*np.exp(-u*u/2)/np.sqrt(2*np.pi)*v
    app = -u*ap*v
    p0 = np.zeros(len(p)); p1 = p0.copy(); p2 = p0.copy()
    for i in range(2):
        c, cp, cpp = a[2*i:2*i+2], ap[2*i:2*i+2], app[2*i:2*i+2]
        s, sp, spp = a[4+i], ap[4+i], app[4+i]
        w, wp, wpp = a[6+i], ap[6+i], app[6+i]
        r = pts-c; r2 = np.sum(r*r, axis=1); rcp = r@cp
        g = w*np.exp(-r2/(2*s*s))
        l1 = wp/w+rcp/s**2+r2*sp/s**3
        l2 = (wpp/w-(wp/w)**2+(-cp@cp+r@cpp)/s**2
              -4*rcp*sp/s**3+r2*spp/s**3-3*r2*sp**2/s**4)
        p0 += g; p1 += g*l1; p2 += g*(l2+l1*l1)
    im = int(np.argmax(p0)); m, mp, mpp = p0[im], p1[im], p2[im]
    return 1.2*prior.V*(p2/m-2*p1*mp/m**2-p0*mpp/m**2+2*p0*mp**2/m**3)


def matmul(a, b):
    COUNTS['matrix_products'] += 1
    return a@b


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cache-root', type=Path, default=ROOT.parent/'b1-finite-feedback/research/finite_feedback_fiber_diagnostic_20261016')
    parser.add_argument('--output', type=Path, default=ROOT/'research/research_f_intrinsic_transport')
    args = parser.parse_args(); start = time.perf_counter()
    out = args.output; out.mkdir(parents=True, exist_ok=True)
    old = args.cache_root; res = old/'results'
    Chart = load_chart_class()
    metrics = json.loads((old/'METRICS.json').read_text())
    hist = {r['anchor']: r for r in metrics['geometry']}
    # Header-only inventory verifies that names such as curvature endpoints do
    # not imply a stored directional Hessian. No paid state arrays are loaded.
    inventory = []
    for p in sorted(res.rglob('*.npz')):
        with zipfile.ZipFile(p) as z:
            keys = [Path(n).stem for n in z.namelist()]
        inventory.append(dict(path=str(p), keys=keys, bytes=p.stat().st_size))
    tests = dict(scope='F0_DECODER_AND_CACHED_J_ONLY', steps=STEPS, rank_rtol=RANK_RTOL,
                 material_metric='[Re chi; Im chi]/24',
                 packing='source-major Re32 then Im32; sqrt(2)/sigma whitening', anchors={},
                 cache_header_inventory=inventory,
                 forbidden_physics_imports=[], new_full_states=0, new_complex_rhs=0,
                 new_lu_factorizations=0, new_mie_calls=0, gpu_training=0)
    rows = []; arrays = {}; report_rows = []
    for key in ('250033_weak', '250034_ring'):
        meta = json.loads((old/'inputs'/f'{key}.json').read_text())
        with np.load(res/f'{key}_CHART.npz') as z: cached = {k: z[k] for k in z.files}
        with np.load(res/f'{key}_JACOBIAN.npz') as z: J = z['J_white_mass']; C = z['C']
        with np.load(old/'inputs'/f'{key}.npz') as z: anchor_chi=z['chi']; sigma=float(z['sigma'])
        COUNTS['npz_files_read'] += 3
        chart = Chart(meta['generator'], key.split('_')[-1]); u = cached['u0']; v=cached['v']; D=cached['D']
        x0, state0 = evaluate(chart, u); DM=mass(D)
        COUNTS['analytic_decoder_first_derivatives'] += 1
        Dcheck=chart.analytic_derivative(u)
        U, s, Vh = np.linalg.svd(DM, full_matrices=False); COUNTS['svds'] += 1
        rank=int(np.sum(s>s[0]*RANK_RTOL)); Q=U[:, :rank]
        A=matmul(J,Q); Ua, sa, _ = np.linalg.svd(A,full_matrices=False); COUNTS['svds'] += 1
        arank=int(np.sum(sa>sa[0]*RANK_RTOL)); Ua=Ua[:, :arank]
        def pout(z): return z-matmul(Ua,matmul(Ua.T,z))
        def pprior(z): return matmul(Q,matmul(Q.T,z))
        w=matmul(DM,v); speed=norm(w)
        d2exact=mass(analytic_second(chart,u,v))
        # Gauge is proportional scaling of the two original Gaussian weights.
        a=chart.lo+(chart.hi-chart.lo)*ndtr(u)
        dadu=(chart.hi-chart.lo)*np.exp(-u*u/2)/np.sqrt(2*np.pi)
        gauge=np.zeros(8); gauge[6:]=a[6:]/dadu[6:]; gauge/=norm(gauge)
        material_gauge_defect=norm(matmul(DM,gauge))/s[0]
        curvature_source=[]
        # Optional explicitly named Hessian cache, never inferred from endpoints.
        for path in (res/f'{key}_CURVATURE.npz', res/f'{key}_HESSIAN.npz'):
            if path.exists():
                with np.load(path) as z:
                    if 'C_Maxwell_white' in z.files: curvature_source=[path,'C_Maxwell_white',z['C_Maxwell_white']]
                    elif 'bmax' in z.files:
                        b=z['bmax']; packed=np.concatenate((b.real.T,b.imag.T),axis=1).ravel()
                        curvature_source=[path,'bmax',packed*np.sqrt(2)/sigma]
                COUNTS['npz_files_read'] += 1
        cm=curvature_source[2] if curvature_source else None
        t=dict(source_chart=str(res/f'{key}_CHART.npz'), source_jacobian=str(res/f'{key}_JACOBIAN.npz'),
            source_input=str(old/'inputs'/f'{key}.npz'), source_generator=str(old/'inputs'/f'{key}.json'),
            anchor_replay_max_abs=float(np.max(abs(x0-anchor_chi))),
            analytic_first_cache_relative_error=rel(Dcheck,D),
            cached_JDg_relative_error=rel(matmul(J,DM),cached['A']),
            prior_coordinate_speed=norm(v), material_speed=speed,
            material_spectrum=s, material_rank=rank, gauge_dimension=8-rank,
            gauge_material_relative_defect=material_gauge_defect,
            gram_effective_condition=float((s[0]/s[rank-1])**2),
            response_material_orthonormal_spectrum=sa, response_rank=arank,
            exact_prior_kernel_dimension=rank-arank, approximate_weak_direction=True,
            exact_fiber_obstruction_instance=False, center=state0,
            maxwell_hessian_status='CACHED' if cm is not None else 'MISSING',
            historical_Maxwell_norm_prior_unit=hist[key]['Maxwell_curvature_norm'],
            historical_decoder_norm_prior_unit=hist[key]['decoder_curvature_norm'],
            hessian_missing_reason='Original executed diagnose.py retained bmax only as local variable; scalar norm persisted in METRICS. No vector/current q key in historical result NPZ headers.',
            finite_difference=[], reparameterization={})
        d2s=[]
        for step in STEPS:
            xp,sp=evaluate(chart,u+step*v); xm,sm=evaluate(chart,u-step*v)
            numerator=xp-2*x0+xm; d2=mass(numerator/step**2); d2s.append(d2)
            ranks=[]
            for uu in (u+step*v,u-step*v):
                ss=np.linalg.svd(mass(chart.analytic_derivative(uu)),compute_uv=False)
                COUNTS['svds'] += 1; COUNTS['analytic_decoder_first_derivatives'] += 1
                ranks.append(int(np.sum(ss>ss[0]*RANK_RTOL)))
            fd=dict(step=step,analytic_D2_relative_error=rel(d2,d2exact),plus=sp,minus=sm,
                side_material_ranks=ranks,
                subtraction_roundoff_scale=float(np.finfo(float).eps*(norm(xp)+2*norm(x0)+norm(xm))/max(norm(numerator),1e-300)),
                max_cell_constant=sp['active_max_cell']==sm['active_max_cell']==state0['active_max_cell'])
            t['finite_difference'].append(fd)
        for index,(label,d2) in enumerate([(f'fd_{h:g}',d) for h,d in zip(STEPS,d2s)]+[('analytic_fixed_max',d2exact)]):
            tangent=pprior(d2); normal=d2-tangent
            cc=matmul(J,tangent); cp=matmul(J,normal); raw=matmul(J,d2)
            oc=pout(cc); op=pout(cp); total=cc+cp
            oldcm=hist[key]['Maxwell_curvature_norm']
            row=dict(anchor=key,derivative=label,prior_speed=norm(v),material_speed=speed,
                material_D2_norm=norm(d2),coordinate_material_D2_norm=norm(tangent),prior_material_D2_norm=norm(normal),
                material_normal_energy_fraction=norm(normal)**2/norm(d2)**2,
                C_coordinate_prior_unit=norm(cc),C_prior_prior_unit=norm(cp),C_raw_prior_unit=norm(raw),
                C_coordinate_material_unit=norm(cc)/speed**2,C_prior_material_unit=norm(cp)/speed**2,
                C_raw_material_unit=norm(raw)/speed**2,
                coordinate_prior_cosine=cosine(cc,cp),
                decoder_cancellation_ratio=norm(total)/(norm(cc)+norm(cp)),
                vector_sum_relative_error=rel(total,raw),
                Qout_coordinate_prior_unit=norm(oc),Qout_prior_prior_unit=norm(op),
                Qout_prior_material_unit=norm(op)/speed**2,
                Qout_coordinate_relative_to_coordinate=norm(oc)/max(norm(cc),1e-300),
                raw_cokernel_equals_prior_relative_error=rel(pout(raw),op),
                historical_Maxwell_norm_prior_unit=oldcm,historical_Maxwell_norm_material_unit=oldcm/speed**2,
                C_Maxwell_status=t['maxwell_hessian_status'],
                C_Maxwell_prior_unit=norm(cm) if cm is not None else 'MISSING',
                C_Maxwell_coordinate_cosine=cosine(cm,cc) if cm is not None else 'MISSING',
                C_Maxwell_prior_cosine=cosine(cm,cp) if cm is not None else 'MISSING',
                total_intrinsic_cokernel_prior_unit=norm(pout(cm+cp)) if cm is not None else 'MISSING',
                full_curvature_vector_reproduction='MISSING' if cm is None else 'DERIVED_NO_OLD_TOTAL_VECTOR',
                historical_norm_only_intrinsic_cokernel_lower=max(0,norm(op)-oldcm),
                historical_norm_only_intrinsic_cokernel_upper=norm(op)+oldcm)
            rows.append(row)
            for name,vec in dict(d2=d2,d2_parallel=tangent,d2_perp=normal,C_coordinate=cc,C_prior=cp,
                C_raw=raw,Qout_C_prior=op,Qout_C_coordinate=oc).items(): arrays[f'{key}__{label}__{name}']=vec
            if label=='analytic_fixed_max': report_rows.append(row)
        t['consecutive_D2_relative_changes']=[rel(d2s[0],d2s[1]),rel(d2s[1],d2s[2])]
        # Same image, invertible at origin: u=u0+z+0.5*a*(v.T z)^2.
        acc=np.eye(8)[0]*.25
        added=matmul(DM,acc)*(v@v)**2
        normal0=d2exact-pprior(d2exact); reexact=d2exact+added
        normal1=reexact-pprior(reexact)
        rp=dict(map='u=u0+z+0.5*(0.25 e0)*(v^T z)^2',
            derivative_at_origin='identity',original_parameter_measure_preserved_for_audit_only=True,
            same_material_velocity_relative_error=0.,
            analytic_intrinsic_material_relative_error=rel(normal1,normal0),
            analytic_intrinsic_cokernel_relative_error=rel(pout(matmul(J,normal1)),pout(matmul(J,normal0))),
            raw_decoder_response_change=norm(matmul(J,added)),finite_difference=[])
        for h in STEPS:
            zp=h*v; zm=-h*v
            up=u+zp+.5*acc*(v@zp)**2; um=u+zm+.5*acc*(v@zm)**2
            xxp,sp=evaluate(chart,up); xxm,sm=evaluate(chart,um)
            red2=mass((xxp-2*x0+xxm)/h**2); ren=red2-pprior(red2)
            rp['finite_difference'].append(dict(step=h,analytic_transformed_D2_relative_error=rel(red2,reexact),
                intrinsic_material_relative_error=rel(ren,normal0),
                intrinsic_cokernel_relative_error=rel(pout(matmul(J,ren)),pout(matmul(J,normal0))),
                active_max_cells=[sp['active_max_cell'],sm['active_max_cell']],
                map_derivative_determinants=[float(1+(v@acc)*(v@zp)),float(1+(v@acc)*(v@zm))]))
        t['reparameterization']=rp
        t['geometry_numerical_checks_pass']=bool(rank==7 and arank==7 and material_gauge_defect<1e-12
            and t['anchor_replay_max_abs']<1e-12 and t['analytic_first_cache_relative_error']<1e-12
            and t['cached_JDg_relative_error']<1e-10
            and all(d['analytic_D2_relative_error']<1e-5 and d['max_cell_constant'] and d['side_material_ranks']==[7,7] for d in t['finite_difference'])
            and rp['analytic_intrinsic_material_relative_error']<1e-10
            and all(d['intrinsic_material_relative_error']<1e-3 for d in rp['finite_difference']))
        arrays[f'{key}__material_tangent_Q']=Q; arrays[f'{key}__material_velocity']=w
        arrays[f'{key}__reparam_D2']=reexact; arrays[f'{key}__reparam_intrinsic_D2']=normal1
        tests['anchors'][key]=t
    # Flat image with arbitrarily curved coordinates: an algebra test, no F solve.
    c=1000.; flatraw=np.array([2*c,0.]); flatQ=np.eye(2)
    flatii=flatraw-flatQ@(flatQ.T@flatraw)
    tests['flat_image_control']=dict(raw_decoder_curvature=norm(flatraw),intrinsic_prior_curvature=norm(flatii),
        geodesic_coordinate_acceleration=[-2*c,0.],corrected_material_acceleration=[0.,0.],pass_test=norm(flatii)==0)
    tests['forbidden_physics_imports']=[m for m in sys.modules if m in ('a30.physics','a25.physics','finite_feedback.diagnose')]
    tests['counts']=dict(COUNTS); tests['wall_seconds']=time.perf_counter()-start
    tests['numerical_geometry_checks_pass']=all(t['geometry_numerical_checks_pass'] for t in tests['anchors'].values()) and not tests['forbidden_physics_imports']
    tests['scientific_gate_decision']='PARENT_REVIEW_REQUIRED'
    with (out/'CURVATURE_DECOMPOSITION.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    arrays['source_note']=np.array('Derived only: original cached CHART/JACOBIAN; decoder-only analytic and fixed-step FD. No exact Maxwell Hessian vector available.')
    np.savez_compressed(out/'GEOMETRY_VECTORS.npz',**arrays)
    (out/'GEOMETRY_TESTS.json').write_text(json.dumps(plain(tests),indent=2,allow_nan=False)+'\n')
    lines=['# F0 原先验材料几何审计','',
        '本次只读取旧 CHART/JACOBIAN 与 generator；decoder-only 求导及矩阵代数。新增 Maxwell states、LU、complex RHS、Mie、GPU 训练均为 **0**。旧 45 states、轨迹、预算账本及 RUN_CLOSED 均未修改。',
        '', '## 度量与定义','',
        '材料坐标为 m(chi)=[Re chi;Im chi]/24，故欧氏范数等于 24×24 材料 RMS。白化数据按每源 Re32/Im32，乘 sqrt(2)/sigma。缓存 J_white_mass 正是这两种度量之间的线性算子。',
        '', '对缓存方向 v（原 Gaussian-CDF 参数 ||v||=1），w=Dg v；SVD 去除共同 Gaussian 权重缩放 gauge 后取材料正交 Q。d2_parallel=QQ^T D²g，d2_perp=(I-QQ^T)D²g。分别计算 C_coordinate=J d2_parallel、C_prior=J d2_perp。二阶量按 ||w||² 归一化，不把参数单位与材料单位混用。',
        '', '这里所称 intrinsic 项是 prior 在固定欧氏材料度量中的外在第二基本形式，不是 Riemann 内蕴曲率张量。去 gauge 不改变原始概率测度，也不创建均匀七维 prior。',
        '', '## 已完成的分解','',
        '| anchor | ||w|| | ||C_coordinate|| | ||C_prior|| | ||raw|| | cosine(coord,prior) | ||Qout C_prior|| |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for r in report_rows:
        lines.append(f"| {r['anchor']} | {r['material_speed']:.9g} | {r['C_coordinate_prior_unit']:.9g} | {r['C_prior_prior_unit']:.9g} | {r['C_raw_prior_unit']:.9g} | {r['coordinate_prior_cosine']:.7g} | {r['Qout_prior_prior_unit']:.9g} |")
    lines += ['', '表中数据响应使用原参数单位；完整逐步和单位材料速度平方结果见 CSV。范数不能加减代替向量分解：coordinate/prior 的 cosine 与 cancellation_ratio 已保存。',
        '', '## 缺失与可用边界','',
        '**两个 anchor 的 exact Maxwell 方向 Hessian 向量均 MISSING。** 原实际执行代码将 bmax/q 留在局部变量，只向 METRICS 写 Maxwell 曲率范数。逐个检查旧 results 全部 NPZ 的 key，未找到该向量或 q；curvature endpoint 只存有限步 state，不能替代方向 Hessian。本次不调用缓存 LU、不解 RHS，也不从粗有限步拟合“exact Hessian”。',
        '', '因而 C_Maxwell 与 C_coordinate/C_prior 的 cosine、三项总和、总 intrinsic 响应和总余核响应均缺失。旧 Maxwell scalar norm 仅提供三角不等式界：max(0,||Qout C_prior||-||C_Maxwell||) ≤ ||Qout(C_Maxwell+C_prior)|| ≤ ||Qout C_prior||+||C_Maxwell||；该界依赖旧 norm 的正确性，不是独立重算。',
        '', '## 差分、gauge 与重参数化','',
        '执行前固定 h=1e-3,5e-4,2.5e-4，全部报告，不筛选最漂亮一步；另独立推导 fixed active-max 的解析 D²g 交叉验证。JSON 给出每步解析误差、相邻变化、减法消去估计、CDF 饱和、active-max cell、侧点 rank。有限点一致不构成整个邻域无边界或恒秩的严格证明。',
        '', '同 image 局部坐标变换 u=u0+z+0.5*(0.25e0)*(v^Tz)^2 在原点导数为 I；固定相同材料速度。raw D²g 增加 Dg(0.25e0)，法向项和其数据余核投影不变。已同时执行解析和三步 FD 对照；映射 Jacobian 行列式在所检查点非零。不声称概率目标在新坐标下自动均匀。',
        '', '两个 prior 材料 rank 都是 7，JQ 也满列秩；当前不存在非零 exact prior tangent kernel。所查缓存 v 是 approximate weak direction，余核二阶响应只能称 fixed-velocity finite-noise 几何诊断，不能称已证实 exact fiber 障碍实例。',
        '', '平坦 image 弯曲坐标控制 g(u,v)=(u+1000v²,v) 给 raw curvature=2000、prior normal curvature=0；geodesic 参数加速度消去材料伪曲率。',
        '', '## 历史边界与交付','',
        'anchor 是公开且已暴露的 TRAIN 已知材料，自己的 GT 就是起点；本次没有独立 anchor、GT 恢复或新反演结果。保留旧 sign-tracking bug 与已发布的轨迹重放说明：此次没有重跑旧轨迹。',
        '', 'geometry numerical checks 与完整科学 gate 分开：本脚本只给数值一致性结果，父线程负责 Gate1 与下一阶段裁决。',
        '', f"新增时间 {tests['wall_seconds']:.6f} 秒；计数 {json.dumps(tests['counts'],ensure_ascii=False)}。",
        '', '来源路径、完整缓存 key inventory、数值误差、条件数和所有 FAIL/MISSING 均在 GEOMETRY_TESTS.json；衍生向量在 GEOMETRY_VECTORS.npz。']
    (out/'INTRINSIC_GEOMETRY_AUDIT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status='NUMERICAL_GEOMETRY_PASS' if tests['numerical_geometry_checks_pass'] else 'CHECK_FAILURE',
        Maxwell_Hessian='MISSING',new_states=0,new_RHS=0,output=str(out),wall_seconds=tests['wall_seconds'])))


if __name__=='__main__': main()
