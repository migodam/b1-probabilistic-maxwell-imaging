"""R12 prediction-only validation/blind rendering, then separate offline scoring.

No shared physical module is imported. The only truth file opening occurs in
score(), after prediction validation and TRUTH_ACCESS_BEGIN.json creation.
Run --help for the small manifest contract. No files are silently overwritten.
"""
import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import label, uniform_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from adaptive_depth import outside_span, orthobasis

HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "STATISTICS_FROZEN.json"
PUBLIC_KEYS = ("materials", "material", "phi", "Phi", "qualified", "noise_qualified", "legal", "layer", "branch", "xs", "params", "x",
               "Fs", "F", "prediction", "prior_cost", "candidate_ids", "initial_weak8", "weaklift",
               "online_selected_index", "occupancy", "occupancies", "occupancyROI", "delta", "deltas", "dmaterial", "ROI", "M")


def clean(obj):
    if isinstance(obj, np.ndarray):
        return clean(obj.tolist())
    if isinstance(obj, np.generic):
        return clean(obj.item())
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, complex):
        return [obj.real, obj.imag]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    return obj


def write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as f:
        json.dump(clean(obj), f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")


def read(path):
    return json.loads(Path(path).read_text())


def csv_write(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows({k: json.dumps(clean(v)) if isinstance(v, (dict, list, np.ndarray)) else clean(v)
                         for k, v in row.items()} for row in rows)


def resolve(path, base):
    p = Path(path)
    return p.resolve() if p.is_absolute() else (base / p).resolve()


def public_npz(path):
    with np.load(path, allow_pickle=False) as z:
        forbidden = [k for k in z.files if ("truth" in k.lower() or k.lower() in
                     ("gt", "xgt", "target_digit", "target_occupancy", "gt_material"))]
        if forbidden:
            raise ValueError("GT-bearing array key in prediction file: " + str(path))
        return {k: np.asarray(z[k]).copy() for k in PUBLIC_KEYS if k in z.files}


def first(z, *names, default=None):
    return next((z[k] for k in names if k in z), default)


def vector_field(value, n, name, dtype=float, default=None):
    if value is None:
        if default is None:
            return None
        return np.full(n, default, dtype=dtype)
    a = np.asarray(value, dtype=dtype)
    if a.ndim == 0:
        a = np.repeat(a, n)
    if a.shape != (n,):
        raise ValueError(name + " must have one value per candidate")
    return a


def load_group(entry, base, cfg):
    g = dict(entry)
    for k in ("dataset", "parent", "method"):
        if k not in g:
            raise ValueError("group requires " + k)
    path = resolve(g["cloud_path"], base) if g.get("cloud_path") else None
    z = public_npz(path) if path else {}
    materials = np.asarray(first(z, "materials", "material", default=np.empty((0, 576))), complex)
    if materials.ndim == 1:
        materials = materials[None]
    if materials.ndim != 2 or materials.shape[1] != 576:
        raise ValueError("full original576 materials required")
    n = len(materials)
    phi = vector_field(first(z, "phi", "Phi"), n, "phi")
    if phi is None and n:
        raise ValueError("candidate data Phi missing")
    phi = np.array([]) if phi is None else phi
    noise_max = float(g.get("noise_phi_max", cfg["noise_Phi_max_default"]))
    finite = np.isfinite(materials).all(axis=1) & np.isfinite(phi)
    legal = vector_field(z.get("legal"), n, "legal", bool)
    claim = vector_field(z.get("qualified"), n, "qualified", bool)
    if claim is None:
        if legal is None:
            raise ValueError("qualified or explicit legal flag required")
        claim = legal & finite & (phi <= noise_max)
    if np.any(claim & (~finite | (phi > noise_max))):
        raise ValueError("claimed qualified candidate fails actual noise/finite gate")
    if legal is not None and np.any(claim & ~legal):
        raise ValueError("claimed qualified but illegal")
    layer = vector_field(z.get("layer"), n, "layer", int, 0)
    xs = first(z, "xs", "params", "x")
    if xs is not None:
        xs = np.asarray(xs, float)
        if xs.ndim == 1:
            xs = xs[None]
        if xs.ndim != 2 or len(xs) != n:
            raise ValueError("candidate parameter matrix shape")
    prior = vector_field(z.get("prior_cost"), n, "prior_cost")
    ids = first(z, "candidate_ids", default=np.arange(n).astype(str))
    if np.asarray(ids).shape != (n,):
        raise ValueError("candidate IDs must align with original immutable order")
    roi = np.asarray(g.get("roi", []), int)
    if g.get("roi_path"):
        rz = public_npz(resolve(g["roi_path"], base))
        roi = np.asarray(rz["ROI"], int)
    anchor = None
    if g.get("anchor_path"):
        az = public_npz(resolve(g["anchor_path"], base))
        am = np.asarray(first(az, "material", "materials"), complex).reshape(576)
        ap = g.get("anchor_phi")
        if ap is None:
            ap = first(az, "phi", "Phi")
        if ap is None and g.get("anchor_status_path"):
            ap = first(read(resolve(g["anchor_status_path"], base)), "phi", "Phi")
        if ap is None:
            raise ValueError("anchor actual data Phi is required")
        anchor = dict(material=am, phi=float(ap), x=first(az, "x", "xs", "params"),
                      occupancy=az.get("occupancy"), delta=az.get("delta"), dmaterial=az.get("dmaterial"),
                      legal=bool(g.get("anchor_legal", True)))
    if anchor is None and n:
        # Only an explicit index can identify an anchor inside a cloud.
        ai = g.get("anchor_index")
        if ai is not None:
            ai = int(ai)
            anchor = dict(material=materials[ai], phi=float(phi[ai]), x=None if xs is None else xs[ai],
                          occupancy=None, delta=None, dmaterial=None, legal=bool(claim[ai]))
    if not roi.size:
        # Known circular ROI is supplied by root; deriving support from a
        # candidate or GT would change the task domain and is forbidden.
        raise ValueError("original fixed ROI must be explicitly provided")
    if roi.shape != (256,) or len(set(roi.tolist())) != 256 or roi.min() < 0 or roi.max() >= 576:
        raise ValueError("original256 ROI indices required")
    occ = first(z, "occupancies", "occupancy", "occupancyROI")
    delta = first(z, "deltas", "delta")
    if occ is not None:
        occ = np.asarray(occ, float)
        if occ.ndim == 1:
            occ = occ[None]
        if occ.shape == (n, 256):
            full = np.zeros((n, 576)); full[:, roi] = occ; occ = full
        if occ.shape != (n, 576):
            raise ValueError("occupancy must be original576 or originalROI256")
    if occ is None and delta is not None:
        ds = vector_field(delta, n, "delta", complex)
        if np.any(ds.real <= 0):
            raise ValueError("positive real increment required for occupancy recovery")
        occ = (materials.real - cfg["host_real"]) / ds.real[:, None]
        if np.max(abs((materials - (cfg["host_real"] + 1j*cfg["host_imag"])) - ds[:, None]*occ), initial=0) > 1e-9:
            raise ValueError("material/delta cannot recover occupancy exactly")
    weak = first(z, "initial_weak8", "weaklift")
    material_reference_basis = None
    if g.get("initial_frame_path"):
        fz = public_npz(resolve(g["initial_frame_path"], base))
        weak = first(fz, "initial_weak8", "weaklift")
        if weak is not None and 'M' in fz:
            material_reference_basis = fz['M'] @ weak
    if weak is not None and (np.asarray(weak).ndim != 2 or np.asarray(weak).shape[1] != 8):
        # A width16 frame cannot be relabeled as initial weak8. Leave this
        # metric unavailable unless the true common reference is supplied.
        g['span_metric_status'] = 'NOT_AVAILABLE: supplied reference is not rank8'
        weak = None
    selected = g.get("online_selected_index", first(z, "online_selected_index"))
    if selected is not None:
        selected = int(np.asarray(selected).item())
        if selected < 0 or selected >= n:
            raise ValueError("sealed online selected index out of range")
    elif claim.any():
        selected = int(np.flatnonzero(claim)[np.argmin(phi[claim])])
    g.update(z=z, materials=materials, phi=phi, qualified=claim, finite=finite, layers=layer, xs=xs,
             prior_cost=prior, ids=np.asarray(ids).astype(str), roi=roi, anchor=anchor,
             occupancy=occ, weak=weak, material_reference_basis=material_reference_basis,
             selected=selected, noise_max=noise_max, source_path=path)
    return g


def load_predictions(manifest_path, config_path=DEFAULT_CONFIG):
    mp = Path(manifest_path).resolve(); m = read(mp); cfg = read(config_path)
    if m.get("predictions_sealed") is not True or m.get("truth_read") is not False:
        raise RuntimeError("separate prediction seal with predictions_sealed=true, truth_read=false required")
    if not cfg.get("frozen_before_R12_truth_scoring"):
        raise RuntimeError("statistics freeze required")
    groups = [load_group(g, mp.parent, cfg) for g in m["groups"]]
    identities = [(g['dataset'], g['parent'], g['method']) for g in groups]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate dataset/parent/method")
    expected = m.get("expected_parents", {})
    for ds in sorted(set(g['dataset'] for g in groups)):
        target = set(expected.get(ds, [g['parent'] for g in groups if g['dataset'] == ds]))
        for method in set(g['method'] for g in groups if g['dataset'] == ds):
            actual = {g['parent'] for g in groups if g['dataset'] == ds and g['method'] == method}
            if actual != target:
                raise ValueError("missing parent/failure group; retain empty groups in manifest")
    return m, groups, cfg


def records(g):
    out = []
    a = g['anchor']
    if a is not None and a['legal'] and np.isfinite(a['material']).all() and a['phi'] <= g['noise_max']:
        out.append(dict(material=a['material'], phi=a['phi'], layer=-1, origin_index=-1, id="ANCHOR", x=a['x']))
    for i in np.flatnonzero(g['qualified']):
        out.append(dict(material=g['materials'][i], phi=g['phi'][i], layer=int(g['layers'][i]),
                        origin_index=int(i), id=str(g['ids'][i]), x=None if g['xs'] is None else g['xs'][i]))
    return out


def unique(records_, tolerance):
    keep = []
    for r in records_:
        if not keep or min(np.linalg.norm(r['material']-k['material'])/24 for k in keep) > tolerance:
            keep.append(r)
    return keep


def topology(mask):
    mask = np.asarray(mask, bool).reshape(24, 24)
    _, components = label(mask, structure=np.array([[0,1,0], [1,1,1], [0,1,0]]))
    inv, n = label(~mask, structure=np.ones((3,3), int))
    edge = set(np.r_[inv[0], inv[-1], inv[:,0], inv[:,-1]].tolist())
    return dict(components=int(components), holes=sum(k not in edge for k in range(1, n+1)))


def ssim(x, y):
    x, y = np.asarray(x).reshape(24,24), np.asarray(y).reshape(24,24)
    ux, uy = uniform_filter(x, 7), uniform_filter(y, 7)
    vx = 49/48*(uniform_filter(x*x,7)-ux*ux)
    vy = 49/48*(uniform_filter(y*y,7)-uy*uy)
    vxy = 49/48*(uniform_filter(x*y,7)-ux*uy)
    m = ((2*ux*uy+.01**2)*(2*vxy+.03**2))/((ux*ux+uy*uy+.01**2)*(vx+vy+.03**2))
    return float(m[3:-3,3:-3].mean())


def structural(pred, truth, threshold):
    p, t = np.asarray(pred).reshape(576), np.asarray(truth).reshape(576)
    pm, tm = p >= threshold, t >= threshold
    inter, union = np.sum(pm & tm), np.sum(pm | tm)
    pt, tt = topology(pm), topology(tm)
    return dict(IoU=float(inter/max(1,union)), Dice=float(2*inter/max(1,pm.sum()+tm.sum())),
                predicted_components=pt['components'], GT_components=tt['components'],
                predicted_holes=pt['holes'], GT_holes=tt['holes'],
                topology_exact=pt == tt, threshold=threshold)


def task_values(materials, roi, cfg):
    x = np.asarray(materials)[:, roi]
    host = cfg['host_real'] + 1j*cfg['host_imag']
    return dict(mean_real=x.real.mean(axis=1), energy=np.mean(abs(x-host)**2,axis=1),
                threshold_count=np.sum(x.real >= cfg['task_pixel_real_threshold'],axis=1))


def cloud_score(rs, truth, g, cfg):
    eps = [cfg['epsilon_primary']] + cfg['epsilon_secondary']
    result = dict(qualified_unique_total=len(rs), qualified_unique_nonanchor=sum(r['origin_index']>=0 for r in rs),
                  available=bool(rs), coverage={str(e):False for e in eps}, nearest_material_MSE=None,
                  nearest_material_RMS=None, nearest_error_right_censored=not bool(rs))
    if not rs:
        return result
    x = np.stack([r['material'] for r in rs]); error = np.mean(abs(x-truth)**2,axis=1)
    nearest = int(np.argmin(error)); phi = np.array([r['phi'] for r in rs])
    dist = np.linalg.norm(x[:,None]-x[None,:],axis=2)/24
    result.update(nearest_material_MSE=float(error[nearest]), nearest_material_RMS=float(np.sqrt(error[nearest])),
                  offline_nearest_member_id=rs[nearest]['id'], coverage={str(e):bool(error[nearest] <= e*e) for e in eps},
                  pair_diameter=float(dist.max()), RMS_dispersion=float(np.sqrt(np.mean(abs(x-x.mean(axis=0))**2))),
                  Phi_min=float(phi.min()), Phi_max=float(phi.max()))
    w = np.exp(-(phi-phi.min())); w /= w.sum()
    result['actions'] = dict(uniform_MSE=float(np.mean(abs(x.mean(axis=0)-truth)**2)),
                             likelihood_weighted_MSE=float(np.mean(abs(w@x-truth)**2)),
                             semantics=cfg['ensemble_semantics'])
    values, target = task_values(x,g['roi'],cfg), task_values(truth[None],g['roi'],cfg)
    result['tasks'] = {}
    thresholds = dict(mean_real=cfg['task_mean_real_event'], energy=cfg['task_energy_event'],
                      threshold_count=cfg['task_threshold_count_event'])
    for name, vals in values.items():
        t = float(target[name][0]); event = vals >= thresholds[name]
        result['tasks'][name] = dict(min=float(vals.min()),max=float(vals.max()),width=float(vals.max()-vals.min()),
                                   GT=t,GT_in_range=bool(vals.min()<=t<=vals.max()),nearest_absolute_error=float(np.min(abs(vals-t))),
                                   event_threshold=thresholds[name],GT_event=t>=thresholds[name],
                                   GT_event_represented=bool(np.any(event==(t>=thresholds[name]))),
                                   uniform_event_score=float(event.mean()),likelihood_event_score=float(w@event))
    return result


def selected_material(g):
    if g['selected'] is not None:
        i=g['selected']
        return g['materials'][i], None if g['occupancy'] is None else g['occupancy'][i], str(g['ids'][i]), bool(g['qualified'][i])
    if g['anchor'] is not None:
        a=g['anchor']; occ=a.get('occupancy')
        if occ is not None and np.asarray(occ).shape==(256,):
            tmp=np.zeros(576);tmp[g['roi']]=occ;occ=tmp
        return a['material'],occ,'ANCHOR_FALLBACK',a['phi']<=g['noise_max'] and a['legal']
    return None,None,'UNAVAILABLE',False


def make_blind(manifest_path, out, config_path=DEFAULT_CONFIG):
    m, groups, cfg = load_predictions(manifest_path,config_path)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(cfg['blind_shuffle_seed']);key={};files=[]
    for di,ds in enumerate(sorted(set(g['dataset'] for g in groups))):
        parents=sorted(set(g['parent'] for g in groups if g['dataset']==ds));parents=list(np.array(parents)[rng.permutation(len(parents))])
        methods=sorted(set(g['method'] for g in groups if g['dataset']==ds));methods=list(np.array(methods)[rng.permutation(len(methods))])
        for mi,method in enumerate(methods):
            code=f'{chr(65+di)}{mi+1:02d}';gg=[next(g for g in groups if (g['dataset'],g['parent'],g['method'])==(ds,p,method)) for p in parents]
            key[code]=dict(dataset=ds,method=method,panels={f'{code}-{j+1:02d}':g['parent'] for j,g in enumerate(gg)})
            for view in ('MATERIAL_REAL','OCCUPANCY'):
                selected=[selected_material(g) for g in gg]
                if view=='OCCUPANCY' and not any(v[1] is not None for v in selected):
                    continue
                nc=min(3,max(1,len(gg)));nr=math.ceil(len(gg)/nc)
                fig,axs=plt.subplots(nr,nc,figsize=(3*nc,3*nr+.6),squeeze=False)
                for j,ax in enumerate(axs.ravel()):
                    ax.axis('off')
                    if j>=len(gg):continue
                    mat,occ,_,qual=selected[j];data=occ if view=='OCCUPANCY' else None if mat is None else mat.real
                    if data is None:ax.text(.5,.5,'UNAVAILABLE',ha='center',va='center')
                    else:ax.imshow(np.asarray(data).reshape(24,24),cmap='gray_r' if view=='OCCUPANCY' else 'viridis',
                                   vmin=0 if view=='OCCUPANCY' else cfg['host_real'],vmax=1 if view=='OCCUPANCY' else 1.9,interpolation='nearest')
                    ax.set_title(f'{code}-{j+1:02d}',fontsize=11)
                fig.suptitle(f'Anonymous {view.lower().replace("_"," ")}; raw 24×24',fontsize=12)
                fig.tight_layout(rect=(0,0,1,.97));path=out/f'BLIND_{code}_{view}.png';fig.savefig(path,dpi=140);plt.close(fig);files.append(str(path.resolve()))
    write(out.parent/'audit/BLIND_KEY_SEPARATE.json',dict(order=key,created_utc=time.time(),GT_read=False,labels_hidden=True,
                                                       note='Do not open until the main reader records the blind visual verdict.'))
    write(out/'BLIND_RENDER_RECEIPT.json',dict(files=files,GT_read=False,new_physics_calls=0,manifest=str(Path(manifest_path).resolve()),
                                            display='nearest-neighbor raw24; common scales; online selection only; no GT-nearest selection'))
    return files


def truth_paths(path):
    p=Path(path).resolve();data=read(p)
    if 'truth' in data:
        return {(r['dataset'],r['parent']):resolve(r['path'],p.parent) for r in data['truth']}
    return {(ds,parent):resolve(v['path'] if isinstance(v,dict) else v,p.parent)
            for ds,parents in data.items() for parent,v in parents.items()}


def score(manifest_path, truth_manifest, out, config_path=DEFAULT_CONFIG, blind_read=None):
    manifest,groups,cfg=load_predictions(manifest_path,config_path)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if any(g['dataset'] in ('N','DIGITS_N') for g in groups):
        if blind_read is None or not Path(blind_read).is_file():
            raise RuntimeError('main-thread blind visual receipt required before digit GT scoring')
        br=read(blind_read)
        if not (br.get('main_thread_viewed') is True or br.get('reader')=='main'):
            raise RuntimeError('blind read receipt must identify actual main-thread viewing')
    paths=truth_paths(truth_manifest)
    if any((g['dataset'],g['parent']) not in paths for g in groups):
        raise ValueError('truth mapping missing expected group')
    write(out/'TRUTH_ACCESS_BEGIN.json',dict(utc=time.time(),manifest=str(Path(manifest_path).resolve()),
          predictions_sealed=True,statistics_id=cfg['id'],role='OFFLINE_ONLY',new_physics_calls=0))
    truths={}
    for identity in dict.fromkeys((g['dataset'],g['parent']) for g in groups):
        # Sole private-truth np.load entry. No model or direction is changed.
        with np.load(paths[identity],allow_pickle=False) as z:
            truth=np.asarray(z['material'],complex).copy()
            occ=np.asarray(z['occupancy'],float).copy() if 'occupancy' in z.files else None
            digit=int(np.asarray(z['target_digit']).item()) if 'target_digit' in z.files else None
        if truth.shape!=(576,) or not np.isfinite(truth).all():raise ValueError('finite original576 GT material required')
        truths[identity]=dict(material=truth,occupancy=occ,digit=digit)
    group_rows=[];candidate_rows=[];layer_rows=[];topology_rows=[];summary={}
    for g in groups:
        identity={k:g[k] for k in ('dataset','parent','method')};gt=truths[(g['dataset'],g['parent'])];rs=unique(records(g),cfg['uniqueness_material_RMS'])
        if gt['digit'] is not None:identity['digit']=gt['digit']
        result=cloud_score(rs,gt['material'],g,cfg);row=dict(**identity,**result)
        row['qualified_unique_RMS002_total'] = len(unique(records(g),cfg['qualification_distinct_RMS']))
        row['noise_qualified_raw_count'] = int(np.sum(g['z'].get('noise_qualified',g['qualified'])))
        selected,occ,sid,squalified=selected_material(g)
        row.update(online_selected_id=sid,online_selected_qualified=squalified,
                   online_selected_MSE=None if selected is None else float(np.mean(abs(selected-gt['material'])**2)),
                   online_selected_RMS=None if selected is None else float(np.sqrt(np.mean(abs(selected-gt['material'])**2))),
                   online_selected_ROI_MSE=None if selected is None else float(np.mean(abs(selected[g['roi']]-gt['material'][g['roi']])**2)))
        if g['anchor'] is not None:row['anchor_MSE']=float(np.mean(abs(g['anchor']['material']-gt['material'])**2))
        group_rows.append(row)
        prev=None
        for lay in sorted(set([-1]+g['layers'].tolist())):
            lr=unique([r for r in records(g) if r['layer']<=lay],cfg['uniqueness_material_RMS'])
            ls=cloud_score(lr,gt['material'],g,cfg);layer_rows.append(dict(**identity,layer=int(lay),**ls))
        for i,(mat,phi,qual,lay,cid) in enumerate(zip(g['materials'],g['phi'],g['qualified'],g['layers'],g['ids'])):
            mse=float(np.mean(abs(mat-gt['material'])**2)) if g['finite'][i] else None
            cr=dict(**identity,index=i,candidate_id=str(cid),layer=int(lay),qualified=bool(qual),Phi=float(phi),
                    prior_cost=None if g['prior_cost'] is None else float(g['prior_cost'][i]),material_MSE=mse,
                    material_RMS=None if mse is None else float(np.sqrt(mse)),
                    anchor_material_RMS=None if g['anchor'] is None else float(np.linalg.norm(mat-g['anchor']['material'])/24))
            if g['weak'] is not None and g['xs'] is not None and g['anchor'] is not None and g['anchor']['x'] is not None:
                cr['parameter_span_outside_initial_weak8']=outside_span(g['xs'][i]-np.asarray(g['anchor']['x']),g['weak'])
            if g['material_reference_basis'] is not None and g['anchor'] is not None:
                dm=mat-g['anchor']['material'];v=np.r_[dm.real,dm.imag]/24
                cr['material_tangent_span_outside_initial_weak8']=outside_span(v,g['material_reference_basis'])
            candidate_rows.append(cr)
        if selected is not None:
            tr=dict(**identity,selected_id=sid,grid='raw24',material_mask=structural(selected.real,gt['material'].real,cfg['material_topology_threshold']))
            if occ is not None and gt['occupancy'] is not None:
                tr.update(occupancy_MSE=float(np.mean((np.asarray(occ)-gt['occupancy'])**2)),
                          ROI_occupancy_MSE=float(np.mean((np.asarray(occ)[g['roi']]-gt['occupancy'][g['roi']])**2)),
                          SSIM_raw24=ssim(occ,gt['occupancy']),occupancy_topology=structural(occ,gt['occupancy'],cfg['occupancy_topology_threshold']))
            else:tr['occupancy_topology_status']='NOT_AVAILABLE; no material amplitude normalization invented'
            topology_rows.append(tr)
    for ds in sorted(set(r['dataset'] for r in group_rows)):
        summary[ds]={}
        for method in sorted(set(r['method'] for r in group_rows if r['dataset']==ds)):
            rr=[r for r in group_rows if r['dataset']==ds and r['method']==method]
            summary[ds][method]=dict(parents=len(rr),available=sum(r['available'] for r in rr),
                coverage_eps075=sum(r['coverage'][str(cfg['epsilon_primary'])] for r in rr),
                coverage_rate=sum(r['coverage'][str(cfg['epsilon_primary'])] for r in rr)/len(rr),
                qualified_unique_mean=float(np.mean([r['qualified_unique_total'] for r in rr])),
                nearest_MSE_mean=None if any(r['nearest_material_MSE'] is None for r in rr) else float(np.mean([r['nearest_material_MSE'] for r in rr])),
                selected_MSE_mean=None if any(r['online_selected_MSE'] is None for r in rr) else float(np.mean([r['online_selected_MSE'] for r in rr])))
    write(out/'GROUP_SCORES.json',group_rows);write(out/'CANDIDATE_SCORES.json',candidate_rows);write(out/'LAYER_SCORES.json',layer_rows);write(out/'TOPOLOGY_SCORES.json',topology_rows)
    csv_write(out/'GROUP_SCORES.csv',[{k:v for k,v in r.items() if k not in ('tasks','actions')} for r in group_rows]);csv_write(out/'CANDIDATE_SCORES.csv',candidate_rows)
    csv_write(out/'TOPOLOGY_SCORES.csv',topology_rows)
    write(out/'SUMMARY.json',dict(status='OFFLINE_FINITE_CLOUD_SCORED',datasets=summary,statistics=cfg,
          new_physics_calls=0,no_post_truth_model_changes=True,historical_P3_visual_status='FAIL',
          scope='Exposed DEV/designed objects as declared by manifest. Nearest GT is offline capacity; no posterior or population guarantee.'))
    plot_scores(group_rows,layer_rows,out)
    write(out/'SCORING_COMPLETE.json',dict(utc=time.time(),new_physics_calls=0,truth_scoring_complete=True))
    return summary


def plot_scores(rows,layers,out):
    for ds in sorted(set(r['dataset'] for r in rows)):
        rr=[r for r in rows if r['dataset']==ds];methods=sorted(set(r['method'] for r in rr));parents=sorted(set(r['parent'] for r in rr))
        fig,axs=plt.subplots(1,3,figsize=(15,4.8));xx=np.arange(len(methods))
        rates=[np.mean([r['coverage']['0.075'] for r in rr if r['method']==m]) for m in methods]
        axs[0].bar(xx,rates,color='#315a88');axs[0].set_ylim(0,1);axs[0].set_ylabel('Fraction of declared parents');axs[0].set_title('ε=.075 material-cloud coverage')
        for pi,p in enumerate(parents):
            vals=[next(r for r in rr if r['parent']==p and r['method']==m) for m in methods]
            axs[1].plot(xx,[r['nearest_material_RMS'] if r['nearest_material_RMS'] is not None else np.nan for r in vals],'.-',alpha=.7,label=p)
            axs[2].plot(xx,[r['qualified_unique_total'] for r in vals],'.-',alpha=.7)
        axs[1].axhline(.075,color='black',linestyle='--',linewidth=1);axs[1].set_title('Offline nearest actual member');axs[1].set_ylabel('Original576 complex material RMS')
        axs[2].set_title('Qualified unique finite-cloud members');axs[2].set_ylabel('Count, anchor included')
        for ax in axs:ax.set_xticks(xx,methods,rotation=45,ha='right');ax.grid(axis='y',alpha=.2)
        if len(parents)<=9:axs[1].legend(fontsize=6,loc='best')
        fig.suptitle(ds+': finite evaluated cloud; exposed DEV/designed objects; no posterior claim',fontsize=11)
        fig.tight_layout();fig.savefig(out/f'{ds}_COVERAGE_WIDTH.png',dpi=150);fig.savefig(out/f'{ds}_COVERAGE_WIDTH.pdf');plt.close(fig)
        fig,axs=plt.subplots(1,2,figsize=(11,4.8))
        for method in methods:
            ls=[r for r in layers if r['dataset']==ds and r['method']==method];ll=sorted(set(r['layer'] for r in ls))
            axs[0].plot(ll,[np.mean([r['qualified_unique_total'] for r in ls if r['layer']==l]) for l in ll],'.-',label=method)
            axs[1].plot(ll,[np.mean([r['coverage']['0.075'] for r in ls if r['layer']==l]) for l in ll],'.-',label=method)
        axs[0].set_ylabel('Mean qualified unique member count');axs[1].set_ylabel('ε=.075 coverage fraction');axs[1].set_ylim(0,1)
        for ax in axs:ax.set_xlabel('Recorded layer prefix (-1 = anchor only)');ax.grid(alpha=.2);ax.legend(fontsize=7)
        fig.suptitle(ds+': layer prefixes use actual candidates; missing/empty groups retained',fontsize=11);fig.tight_layout();fig.savefig(out/f'{ds}_LAYER_PREFIX.png',dpi=150);plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['validate','blind','score']);ap.add_argument('--manifest',required=True);ap.add_argument('--out',required=True)
    ap.add_argument('--config',default=str(DEFAULT_CONFIG));ap.add_argument('--truth-manifest');ap.add_argument('--blind-read');a=ap.parse_args()
    if a.action=='validate':
        _,gg,_=load_predictions(a.manifest,a.config);print(json.dumps(dict(status='PREDICTIONS_VALIDATED_NO_TRUTH_READ',groups=len(gg),new_physics_calls=0)));return
    if a.action=='blind':print(json.dumps(make_blind(a.manifest,a.out,a.config)));return
    if not a.truth_manifest:raise RuntimeError('separate --truth-manifest required only for score')
    print(json.dumps(score(a.manifest,a.truth_manifest,a.out,a.config,a.blind_read),indent=2))


if __name__=='__main__':
    main()
