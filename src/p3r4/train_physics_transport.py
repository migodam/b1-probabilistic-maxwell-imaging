"""Released full66 explicit-density refinement with actual Maxwell adjoint VJPs.

No paid action is performed on import. This entry point owns a physics Book
and both campaign physics/training locks; never nest it in run_r4_job.py.
Only a separately frozen, explicitly released JSON permits execution. The
pure-supervised arm pays the same fixed physical validation roster but no
training physics. Existing inference checkpoint formats remain unchanged.
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np
import torch

from .anchor_density import context_features, base_parameters, FORMAT as ANCHOR_FORMAT
from .conditional_coordinates import ConditionalCoordinateMap
from .core import OUT, Book, BudgetStop, Target, from_probit, probit_derivative, pack, save, write
from .density import make_model, normal_log_prob
from .train_density import load_data
from .train_anchor_density import load_anchor_data

RAW_FORMAT = 'p3r4_full66_density_v1'
RELEASE_FORMAT = 'p3r4_physics_transport_release_v1'


def identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def verify_identity(record):
    if identity(record['path']) != record:
        raise ValueError('frozen file byte/mtime identity changed: '+record['path'])


def unpack_observations(packed):
    """Inverse of core.pack: per source32 real then32 imaginary samples."""
    packed = np.asarray(packed)
    if packed.ndim != 2 or packed.shape[1] != 512 or np.iscomplexobj(packed) or not np.isfinite(packed).all():
        raise ValueError('finite real source-major y[B,512] required')
    source = packed.reshape(-1, 8, 64)
    return (source[:, :, :32]+1j*source[:, :, 32:]).transpose(0, 2, 1)


class ActualPhi(torch.autograd.Function):
    """Forward pays native Full+adjoint; backward uses only saved CPU VJP.

    The oracle returns dPhi/dx66, including both probit property derivatives.
    No finite-difference, surrogate, decoder-only, or stale-anchor gradient is
    substituted. First derivatives are supported; double backward is not.
    """
    @staticmethod
    def forward(ctx, x, observations, oracle):
        if x.ndim != 2 or x.shape[1] != 66 or len(x) > 4:
            raise ValueError('actual physics microbatch must have1..4 full66 points')
        values, gradients = oracle.value_and_gradient(x.detach().cpu().double().numpy(), observations)
        values, gradients = np.asarray(values), np.asarray(gradients)
        if values.shape != (len(x),) or gradients.shape != (len(x), 66) or not np.isfinite(values).all() or not np.isfinite(gradients).all():
            raise FloatingPointError('nonfinite/malformed physical likelihood VJP')
        # Retain only the small detached gradient, never a Point, LU or graph.
        ctx.save_for_backward(torch.as_tensor(gradients.copy(), dtype=x.dtype, device=x.device))
        return torch.as_tensor(values.copy(), dtype=x.dtype, device=x.device)

    @staticmethod
    @torch.autograd.function.once_differentiable
    def backward(ctx, grad_output):
        gradient, = ctx.saved_tensors
        return grad_output[:, None]*gradient, None, None


class TargetLikelihoodOracle:
    """Serial native Target calls with immediate Point release and8+8 RHS.

    Target itself persists paid endpoints/derivatives. This object stores no
    Points and clears only receiver-transfer cache references between queries.
    """
    def __init__(self, target, costs, max_rhs, query_directory=None):
        self.target, self.costs, self.max_rhs = target, costs, int(max_rhs)
        self.query_directory = Path(query_directory) if query_directory is not None else target.book.jobdir/'transport_queries'
        self.query_directory.mkdir(parents=True,exist_ok=True)
        self.query_count = 0

    def _one(self, x, y, gradient, row):
        x,y = np.asarray(x),np.asarray(y)
        if x.shape != (66,) or not np.isfinite(x).all() or y.shape != (32,8) or y.dtype != np.complex128 or not np.isfinite(y).all():
            raise ValueError('finite raw x66 and ORIGINAL complex128 observation required')
        required = 16 if gradient else 8
        if self.costs['RHS']+required > self.max_rhs:
            raise BudgetStop('physics transport cannot start an incomplete forward/VJP pair')
        before = self.costs['RHS']
        self.query_count += 1
        folder = self.query_directory/f'q{self.query_count:08d}'
        folder.mkdir()
        context = {k:v for k,v in self.target.book.context.items() if k not in ('endpoint','purpose')}
        context['query_row_in_batch'] = row
        if 'parents' in context:
            context['parent'] = str(context['parents'][row])
        if 'sample_indices' in context:
            context['sample_index'] = int(context['sample_indices'][row])
        # This happens BEFORE Target.evaluate. Native q alone cannot recover
        # finite eta when a floating-point CDF rounds to an endpoint.
        save(folder/'BEFORE.npz',x=x.copy(),observation=y.copy(),q=from_probit(x),
             parent=np.array(context.get('parent','')),query_index=np.array(self.query_count))
        write(folder/'BEFORE.json',dict(status='SAVED_BEFORE_FULL',context=context,
              gradient_requested=gradient,expected_RHS=required,RHS_before=before,
              observation='original complex128; frozen source provenance in fit contract'))
        point = None
        value = derivative = None
        query_status,query_error = 'COMPLETE',None
        try:
            point = self.target.evaluate(from_probit(x), purpose='transport_train_VJP' if gradient else 'transport_fixed_validation')
            value = self.target.phi(point, y)
            derivative = np.asarray(self.target.phi_gradient(point, y))*probit_derivative(x) if gradient else None
            if self.costs['RHS']-before != required:
                raise AssertionError('actual forward/adjoint RHS differs from frozen8/16 contract')
            if not np.isfinite(value) or (gradient and not np.isfinite(derivative).all()):
                raise FloatingPointError('nonfinite actual physical likelihood or gradient')
            return value, derivative
        except BaseException as exc:
            query_status,query_error = 'FAILED',repr(exc)
            raise
        finally:
            try:
                endpoint = point.key if point is not None else self.target.book.context.get('endpoint','')
                save(folder/'AFTER.npz',Phi=np.array(value if value is not None else np.nan),
                     gradient_x=np.asarray(derivative) if derivative is not None else np.empty(0),
                     endpoint=np.array(endpoint),paid_RHS=np.array(self.costs['RHS']-before))
                write(folder/'AFTER.json',dict(status=query_status,error=query_error,context=context,
                      endpoint=endpoint,paid_RHS=self.costs['RHS']-before,gradient_available=derivative is not None))
            finally:
                del point
                if self.target.model is not None:
                    self.target.model._transfers.clear()

    def value_and_gradient(self, x, observations):
        x, observations = np.asarray(x), np.asarray(observations)
        if x.shape != (len(observations), 66) or observations.shape[1:] != (32, 8) or not 1 <= len(x) <= 4:
            raise ValueError('aligned bounded x[B,66], complex observations[B,32,8] required')
        values, gradients = [], []
        for row,(a, y) in enumerate(zip(x, observations)):
            value, derivative = self._one(a, y, True,row)
            values.append(value); gradients.append(derivative)
        return np.asarray(values), np.asarray(gradients)

    def values(self, x, observations):
        if len(x) > 4:
            raise ValueError('bounded physical validation microbatch<=4 required')
        if np.asarray(x).shape != (len(observations),66) or np.asarray(observations).shape[1:] != (32,8):
            raise ValueError('aligned raw x and original observation roster required')
        return np.asarray([self._one(a, y, False,row)[0] for row,(a, y) in enumerate(zip(x, observations))])


class DifferentiableDensity:
    """Existing explicit model + differentiable frozen raw-coordinate maps.

    All map constants depend only on y-derived context and TRAIN normalizers.
    Sampling returns its actual inverse-map log-Jacobian without detaching x.
    Both direct and inverse density paths are with respect to raw dx66.
    """
    def __init__(self, payload, device='cpu', dtype=torch.float32):
        if payload.get('format') not in (RAW_FORMAT, ANCHOR_FORMAT) or payload.get('kind') not in ('spline', 'diag_gaussian'):
            raise ValueError('only pretrained explicit full66 spline/diagonal Gaussian checkpoints are released')
        if payload.get('model_config',{}).get('dimension') != 66:
            raise ValueError('all66 stochastic coordinates must be explicit')
        self.payload, self.stats = copy.deepcopy(payload), payload['normalization']
        self.device, self.dtype = torch.device(device), dtype
        self.model = make_model(payload['kind'], **payload['model_config']).to(device=self.device, dtype=dtype)
        self.model.load_state_dict(payload['state_dict'])
        self.anchor = payload['anchor_contract'] if payload['format'] == ANCHOR_FORMAT else None

    def tensor(self, value):
        return torch.as_tensor(np.asarray(value), dtype=self.dtype, device=self.device)

    def prepare(self, data, ids):
        ids = np.asarray(ids)
        y = data['y'][ids]
        count = len(y)
        if self.anchor is None:
            features, center = y, np.zeros((count, 66))
            matrix = inverse_rows = np.broadcast_to(np.eye(66), (count, 66, 66))
            ld = np.zeros(count)
        else:
            fields = {k: data[k][ids] if data.get(k) is not None else None for k in ('context_x', 'context_J', 'base_mean', 'base_chol')}
            ac = self.anchor
            features = context_features(y, fields['context_x'], fields['context_J'], ac['j_summary'])
            mean, _ = base_parameters(fields['context_x'], fields['base_mean'], fields['base_chol'], ac['anchor_scale'])
            mode = ac.get('coordinate_preconditioner','none')
            if mode != 'none' and not ac['residual_center']:
                raise ValueError('GN/random map requires anchor residual centering')
            center = fields['context_x'] if mode != 'none' else (mean if ac['residual_center'] else np.zeros_like(mean))
            coordinate = ConditionalCoordinateMap(mode, fields['context_x'], fields['context_J'], ac.get('alignment_rotation'))
            matrix, ld = coordinate.matrix(), coordinate.logabsdet
            inverse_rows = coordinate.inverse(np.broadcast_to(np.eye(66), (count, 66, 66)))
        stats = self.stats
        return dict(y=self.tensor((features-stats['y_mean'])/stats['y_scale']), center=self.tensor(center),
                    matrix=self.tensor(matrix.copy()), inverse_rows=self.tensor(inverse_rows.copy()),
                    logdet=self.tensor(ld-np.log(stats['x_scale']).sum()+stats.get('frame_logabsdet', 0.)),
                    frame=self.tensor(stats['frame']), frame_inverse=self.tensor(stats.get('frame_inverse', stats['frame'].T)),
                    mean=self.tensor(stats['x_mean']), scale=self.tensor(stats['x_scale']))

    def log_prob(self, x, prepared):
        delta = torch.einsum('bij,bj->bi', prepared['matrix'], x-prepared['center'])
        standardized = (delta@prepared['frame']-prepared['mean'])/prepared['scale']
        return self.model.log_prob(standardized, prepared['y'])+prepared['logdet']

    def reparameterize(self, eps, prepared):
        standardized, network_logdet = self.model.transform(eps, prepared['y'])
        residual = (standardized*prepared['scale']+prepared['mean'])@prepared['frame_inverse']
        x = torch.einsum('bi,bij->bj', residual, prepared['inverse_rows'])+prepared['center']
        logq = normal_log_prob(eps)-network_logdet+prepared['logdet']
        return x, logq


def lambda_at(update, schedule):
    for row in schedule:
        if row['start'] <= update <= row['end']:
            return float(row['value'])
    raise ValueError('lambda schedule does not cover update')


def validation_updates(updates, cadence):
    return sorted({0, updates, *range(cadence, updates+1, cadence)})


def observation_provenance(directory, parents):
    """Metadata only: freeze allowlisted original files without reading GT."""
    directory = Path(directory).resolve()
    return dict(identity=identity(directory),role='original_observation_only',
        schema=['observation','target_id','sigma'],dtype='complex128',shape=[32,8],
        files={str(parent):identity(directory/(str(parent)+'.npz')) for parent in parents})


def verify_observation_provenance(record, parents, *, canonical=True):
    directory = Path(record['identity']['path']).resolve()
    if canonical and directory != (OUT/'inputs/online').resolve():
        raise ValueError('physical likelihood must use original OUT/inputs/online observations')
    verify_identity(record['identity'])
    if record.get('role') != 'original_observation_only' or record.get('schema') != ['observation','target_id','sigma'] or record.get('dtype') != 'complex128' or record.get('shape') != [32,8]:
        raise ValueError('strict original observation provenance contract required')
    if set(record['files']) != set(parents):
        raise ValueError('original observation identity roster differs from frozen1056 parents')
    for parent in parents:
        source = record['files'][parent]
        if Path(source['path']).resolve() != directory/(parent+'.npz'):
            raise ValueError('original observation filename/parent mismatch')
        verify_identity(source)


def load_original_observations(release, train, val, train_ids, val_ids):
    """Cache only1056 original Y arrays; no raw-field/current/LU/GT access."""
    record = release['observation_directory']
    parents = release['physics_train_parents']+release['validation_parents']
    # Canonical-directory identity was already checked by validate_release;
    # check the exact source files again around their actual reads.
    verify_observation_provenance(record,parents)
    target = release['target_config']
    cached = {}
    for data,indices in ((train,train_ids),(val,val_ids)):
        for index in indices:
            parent = str(data['parent'][index])
            source = record['files'][parent]
            with np.load(source['path'],allow_pickle=False) as z:
                if set(z.files) != {'observation','target_id','sigma'}:
                    raise ValueError('observation file contains non-allowlisted arrays: '+parent)
                y = z['observation']
                if y.shape != (32,8) or y.dtype != np.complex128 or not np.isfinite(y).all():
                    raise ValueError('original complex128 observation required: '+parent)
                if z['target_id'].shape != () or str(z['target_id']) != target['target_id'] or z['sigma'].shape != () or float(z['sigma']) != target['noise']['sigma']:
                    raise ValueError('original observation target/sigma mismatch: '+parent)
                if not np.array_equal(pack(y).astype(np.float32),np.asarray(data['y'][index],np.float32)):
                    raise ValueError('NN y source-major float32 does not match original observation: '+parent)
                cached[parent] = y.copy()
            verify_identity(source)
    verify_identity(record['identity'])
    if set(cached) != set(parents):
        raise ValueError('incomplete original observation roster; no rounded-y fallback')
    return cached


def original_batch(cached, parents):
    return np.stack([cached[str(parent)] for parent in parents])


def persist_transport_batch(folder, phase, raw_x, logq, prior, parents, update, lam,
                            query_indices, *, phi=None, error=None):
    """Persist actual sampled coordinates/density before any physical call."""
    folder=Path(folder)
    save(folder/(phase+'.npz'),x=np.asarray(raw_x),logq_x=np.asarray(logq),
         prior_energy=np.asarray(prior),parent=np.asarray(parents).astype(str),
         update=np.array(update),lambda_value=np.array(lam),query_index=np.asarray(query_indices),
         Phi=np.asarray(phi) if phi is not None else np.full(len(raw_x),np.nan))
    write(folder/(phase+'.json'),dict(status='FAILED' if error is not None else
          ('SAVED_BEFORE_FULL' if phase=='BEFORE' else 'COMPLETE'),error=error,
          update=update,lambda_value=lam,reference_measure='raw Gaussianized dx66',
          no_candidate_selection=True,query_indices=np.asarray(query_indices).tolist()))


def release_template(checkpoint, train, validation, seeds):
    """Create reviewable inert metadata; caller must explicitly release it."""
    train_parents=[f'R2_train_{i:05d}' for i in range(1024)]
    val_parents=[f'R2_validation_{i:05d}' for i in range(32)]
    return dict(format=RELEASE_FORMAT, released=False,
        objective='inclusive_NLL_plus_lambda_reverse_ELBO_same_full66_target',
        checkpoint=identity(checkpoint), train=identity(train), validation=identity(validation),
        target_config=json.loads((OUT/'configs/TARGET.json').read_text()), seeds=list(seeds), devices=['cpu','mps'],
        physics_train_parents=train_parents,validation_parents=val_parents,
        observation_directory=observation_provenance(OUT/'inputs/online',train_parents+val_parents),
        array_codec='stored',inclusive_NLL_batch=128,
        updates=3000, microbatch=4, validation_samples=4, validation_noise_seed=20261025142,
        eval_every=300, validation_lambda=.1, lr=5e-5, grad_clip=5, weight_decay=0.,
        lambda_schedule=[dict(start=1,end=300,value=.001),dict(start=301,end=900,value=.01),
                         dict(start=901,end=1800,value=.1),dict(start=1801,end=3000,value=1.)],
        max_rhs=dict(physics=256000,supervised=20000), max_wall=3600, rss_limit_mib=2048)


def validate_release(release, args):
    if release.get('format') != RELEASE_FORMAT or release.get('released') is not True:
        raise ValueError('explicit frozen/released physics-transport contract required')
    if release.get('objective') != 'inclusive_NLL_plus_lambda_reverse_ELBO_same_full66_target':
        raise ValueError('objective declaration mismatch')
    if args.seed not in release['seeds'] or args.mode not in ('physics', 'supervised'):
        raise ValueError('seed/mode not released')
    for key in ('checkpoint', 'train', 'validation'):
        verify_identity(release[key])
        if any(word in release[key]['path'].lower() for word in ('locked', 'offline_only')):
            raise ValueError('locked/offline-only inputs forbidden')
    if release['microbatch'] != 4 or release['validation_samples'] != 4 or release['validation_noise_seed'] != 20261025142:
        raise ValueError('frozen microbatch4/K4/common validation seed required')
    if release['validation_parents'] != [f'R2_validation_{i:05d}' for i in range(32)]:
        raise ValueError('fixed first32 independent VAL parents required')
    if release['physics_train_parents'] != [f'R2_train_{i:05d}' for i in range(1024)]:
        raise ValueError('physics TRAIN pool must be frozen first1024; GT endpoints never define proposals')
    verify_observation_provenance(release['observation_directory'],release['physics_train_parents']+release['validation_parents'])
    if release.get('inclusive_NLL_batch') != 128 or release.get('array_codec') != 'stored':
        raise ValueError('frozen inclusive NLL batch128 and per-process stored array codec required')
    if release['updates'] != 3000 or release['eval_every'] != 300 or release['validation_lambda'] != .1:
        raise ValueError('released3000/every300/fixed composite lambda.1 protocol required')
    required_schedule = [dict(start=1,end=300,value=.001),dict(start=301,end=900,value=.01),dict(start=901,end=1800,value=.1),dict(start=1801,end=3000,value=1.)]
    if release['lambda_schedule'] != required_schedule or release['lr'] != 5e-5 or release['grad_clip'] != 5:
        raise ValueError('frozen lambda/LR/gradient-clip mismatch')
    event_count = len(validation_updates(release['updates'], release['eval_every']))
    validation_rhs = event_count*32*4*8
    training_rhs = release['updates']*4*16 if args.mode == 'physics' else 0
    cap = release['max_rhs'][args.mode]
    if cap > (256000 if args.mode == 'physics' else 20000) or cap < validation_rhs+training_rhs:
        raise ValueError('physical training/validation roster does not fit released hard RHS cap')
    if release['max_wall'] <= 15 or release.get('rss_limit_mib', 2048) > 3072:
        raise ValueError('bounded positive wall and<=3GiB RSS release required')
    return dict(training_RHS=training_rhs, validation_RHS=validation_rhs, planned_RHS=training_rhs+validation_rhs,
                validation_updates=validation_updates(release['updates'], release['eval_every']))


def load_inputs(payload, release):
    if payload['format'] == ANCHOR_FORMAT:
        ac = payload['anchor_contract']
        loader = lambda path: load_anchor_data(path, ac['j_summary'], 'none', ac.get('coordinate_preconditioner', 'none'))
    else:
        loader = load_data
    train, val = loader(release['train']['path']), loader(release['validation']['path'])
    if train['parent'] is None or val['parent'] is None or set(train['parent']) & set(val['parent']):
        raise ValueError('explicit TRAIN/VAL disjoint independent parent IDs required')
    expected = (1024,256) if payload['format'] == ANCHOR_FORMAT else (32768,2048)
    if (len(train['x']), len(val['x'])) != expected:
        raise ValueError('raw32768/2048 or anchor1024/256 dataset contract required')
    # The pretrained normalizers stay frozen. No VAL/TRAIN refitting occurs.
    val_index = {parent:i for i,parent in enumerate(val['parent'])}
    train_index = {parent:i for i,parent in enumerate(train['parent'])}
    try:
        vi = np.asarray([val_index[p] for p in release['validation_parents']])
        ti = np.asarray([train_index[p] for p in release['physics_train_parents']])
    except KeyError as exc:
        raise ValueError('released parent absent from declared data panel') from exc
    return train, val, ti, vi


def run(args):
    release = json.loads(args.release.read_text())
    plan = validate_release(release, args)
    if args.device not in release['devices']:
        raise ValueError('device not released')
    if args.device == 'mps' and not torch.backends.mps.is_available():
        raise ValueError('MPS requested but unavailable')
    if os.environ.get('R4_JOBDIR'):
        raise ValueError('do not nest a physical training Book inside analysis-only runner')
    if args.out.exists():
        raise FileExistsError('fresh output required; no automatic resume/replay')
    if not args.job or '/' in args.job:
        raise ValueError('explicit unique job name required')
    # Only this standalone job process changes codec; the release fixes it
    # identically for VI and supervised controls. Book records the value.
    os.environ['R4_ARRAY_CODEC'] = release['array_codec']
    torch.set_num_threads(1)
    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    noise = torch.Generator().manual_seed(args.seed+1103)
    # Acquire the extra training lock first; Book acquires physics.lock.
    extra_lock = (OUT/'training.lock').open('a')
    book = None
    try:
        fcntl.flock(extra_lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        book = Book('I', args.job, max_wall=release['max_wall'], gpu=args.device=='mps', kind='physics')
        _run_released(args, release, plan, book, rng, noise)
    finally:
        if book is not None and not book.finished:
            book.close('FAILED', 'unexpected exit before terminal training status')
        fcntl.flock(extra_lock, fcntl.LOCK_UN)
        extra_lock.close()


def _run_released(args, release, plan, book, rng, noise):
    started, cpu_started = time.monotonic(), time.process_time()
    reverse_rng = np.random.default_rng(args.seed+2237)
    costs = dict(RHS=0, states=0, LU=0)
    original_update = book._update
    max_rhs = release['max_rhs'][args.mode]
    def limited(event, **kw):
        if costs['RHS']+kw.get('rhs', 0) > max_rhs:
            raise BudgetStop('physics transport fit hard RHS cap')
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        # Always permit receipt closure even after a memory stop.
        if event != 'JOB_END' and rss > release.get('rss_limit_mib', 2048)*1024**2:
            raise BudgetStop('physics transport fit RSS cap')
        result = original_update(event, **kw)
        for k, source in (('RHS','rhs'),('states','states'),('LU','lu')):
            costs[k] += kw.get(source, 0)
        return result
    book._update = limited
    status, error, history = 'COMPLETE', None, []
    contract = {}
    chosen_state, chosen_value, chosen_update = None, float('inf'), None
    completed, train_metrics = 0, None
    args.out.mkdir(parents=True)
    try:
        write(book.jobdir/'COMMAND.json', dict(argv=sys.argv, cwd=str(Path.cwd()), physical_calls_allowed=True,
              GPU_allocation='whole job wall' if args.device=='mps' else 'none', locks=['physics.lock','training.lock']))
        write(args.out/'FROZEN_RELEASE.json', release)
        payload = torch.load(release['checkpoint']['path'], map_location='cpu', weights_only=False)
        density = DifferentiableDensity(payload, args.device)
        train, val, physics_ids, validation_ids = load_inputs(payload, release)
        for key in ('checkpoint','train','validation'):
            verify_identity(release[key])
        observations = load_original_observations(release,train,val,physics_ids,validation_ids)
        optimizer = torch.optim.AdamW(density.model.parameters(), lr=release['lr'], weight_decay=release.get('weight_decay', 0.))
        target = Target(book, keep_fields=True)
        if target.cfg != release['target_config']:
            raise ValueError('Full statistical target changed since release')
        oracle = TargetLikelihoodOracle(target, costs, max_rhs)
        fixed_eps = torch.randn((128,66), generator=torch.Generator().manual_seed(release['validation_noise_seed']))
        # Fixed order: parent0 samples0..3, then parent1, etc.
        validation_repeated = np.repeat(validation_ids, 4)
        contract = dict(format=payload['format'], kind=payload['kind'], seed=args.seed, mode=args.mode,
            model_config=payload['model_config'], parameter_count=sum(p.numel() for p in density.model.parameters()),
            train_size=len(train['x']), validation_size=len(val['x']), physics_train_size=len(physics_ids),
            updates=release['updates'], microbatch=4, inclusive_NLL_batch=128,
            learning_rate=release['lr'], gradient_clip=release['grad_clip'],
            lambda_schedule=release['lambda_schedule'], validation_lambda=.1, criterion='validation_full_NLL_plus_point1_reverse_ELBO',
            source_checkpoint=release['checkpoint'], preserved_normalization=True, reference_measure='raw Gaussianized dx66',
            physical_queries='actual unchanged Full TM24,8 source forward RHS,8 source adjoint RHS per training sample',
            original_observation_provenance=release['observation_directory'],
            physical_observation='strict allowlisted original complex128 Y; rounded NN y is never used by Full likelihood',
            array_codec=release['array_codec'],persist_current_field=True,
            query_persistence='rawx66+originalY before every Full; Phi/gradient_x/endpoint/failure afterward; perbatch current logq and fixed-val raw draws retained',
            train_physics='NONE' if args.mode=='supervised' else 'samples from current exact q given released TRAIN y/context only',
            validation='all proper NLL; fixed32 VAL parents K4 fixed base noise; actual8-RHS Full per draw; same roster for both arms',
            control_cost='pure supervised control pays validation only; actual cost Pareto comparison, no RHS padding',
            objective='inclusive NLL + scheduled lambda*(Phi+.5||x||²+logq_x); same target; omitted constants fixed in theta',
            deployment='existing measurement-only/explicit-y-derived-anchor predictor; no GT, corrected-sample density or resampling wrapper',
            finite_fit_posterior_qualification='NOT_ESTABLISHED', update_zero_eligible=True, plan=plan,
            max_rhs=max_rhs, device=args.device, book_receipt=str(book.jobdir/'CLOSED.json'))
        write(args.out/'RUN_CONTRACT.json', contract)
        save(book.jobdir/'transport_batches/FIXED_VALIDATION_BASE.npz',epsilon=fixed_eps.numpy(),
             parent=np.repeat(val['parent'][validation_ids],4),sample_index=np.tile(np.arange(4),32))

        def checkpoint(state, update, criterion):
            result = copy.deepcopy(payload)
            result.update(state_dict=state, seed=args.seed, selected_update=update, selection_metric=criterion,
                          criterion=contract['criterion'], contract={**payload.get('contract',{}), **contract},
                          refinement_contract=contract, posterior_fidelity='NOT_ESTABLISHED')
            return result

        def evaluate(update):
            nonlocal chosen_state, chosen_value, chosen_update
            density.model.eval()
            nll_sum = 0.
            with torch.no_grad():
                for start in range(0, len(val['x']), 16):
                    ii = np.arange(start,min(start+16,len(val['x'])))
                    prep = density.prepare(val, ii)
                    # MPS has no float64 dtype. Transfer the same computed
                    # float32 log densities first, then reduce on CPU float64.
                    nll_sum += float((-density.log_prob(density.tensor(val['x'][ii]), prep)).cpu().double().sum())
                components,raw_draws = [],[]
                for start in range(0,128,4):
                    ii = validation_repeated[start:start+4]
                    prep = density.prepare(val, ii)
                    x, logq = density.reparameterize(fixed_eps[start:start+4].to(density.device), prep)
                    raw = x.cpu().double().numpy()
                    book.context = dict(update=update, role='fixed_validation', parents=val['parent'][ii].tolist(), sample_indices=list(range(4)))
                    prior = .5*np.square(raw).sum(1)
                    lq = logq.cpu().double().numpy()
                    batch_folder=book.jobdir/f'transport_batches/validation/u{update:05d}/b{start//4:03d}'
                    queries=np.arange(oracle.query_count+1,oracle.query_count+5)
                    persist_transport_batch(batch_folder,'BEFORE',raw,lq,prior,val['parent'][ii],update,.1,queries)
                    try:
                        phi = oracle.values(raw, original_batch(observations,val['parent'][ii]))
                    except BaseException as exc:
                        persist_transport_batch(batch_folder,'AFTER',raw,lq,prior,val['parent'][ii],update,.1,queries,error=repr(exc))
                        raise
                    persist_transport_batch(batch_folder,'AFTER',raw,lq,prior,val['parent'][ii],update,.1,queries,phi=phi)
                    components.append(np.column_stack((phi,prior,lq)))
                    raw_draws.append(raw)
            components = np.concatenate(components)
            raw_draws = np.concatenate(raw_draws)
            nll = nll_sum/len(val['x'])
            reverse = float(components.sum(1).mean())
            criterion = nll+.1*reverse
            if not np.isfinite(criterion):
                raise FloatingPointError('nonfinite composite validation; do not drop fixed failed slots')
            row = dict(update=update, validation_full_NLL=nll, validation_Phi=float(components[:,0].mean()),
                       validation_prior_energy=float(components[:,1].mean()), validation_logq=float(components[:,2].mean()),
                       validation_entropy=float(-components[:,2].mean()), validation_reverse_ELBO=reverse,
                       validation_composite=criterion, train=train_metrics, costs=costs.copy(),
                       wall_seconds=time.monotonic()-started, cpu_seconds=time.process_time()-cpu_started)
            save(args.out/f'validation_components/u{update:05d}.npz',parent=val['parent'][validation_ids],
                 x=raw_draws.reshape(32,4,66),logq_x=components[:,2].reshape(32,4),
                 components=components.reshape(32,4,3),component_names=np.array(['Phi','prior_energy','logq']),full_NLL=np.array(nll))
            history.append(row)
            with (args.out/'loss_curve.jsonl').open('a') as stream:
                stream.write(json.dumps(row,allow_nan=False)+'\n')
            state = {k:v.detach().cpu().clone() for k,v in density.model.state_dict().items()}
            if criterion < chosen_value:
                chosen_value, chosen_update, chosen_state = criterion, update, state
                torch.save(checkpoint(state, update, criterion), args.out/'BEST.pt')
                save(args.out/'selected_validation_components.npz', parent=val['parent'][validation_ids],
                    x=raw_draws.reshape(32,4,66),logq_x=components[:,2].reshape(32,4),
                    components=components.reshape(32,4,3), component_names=np.array(['Phi','prior_energy','logq']), full_NLL=np.array(nll))
            torch.save(dict(**checkpoint(state, update, criterion), optimizer=optimizer.state_dict(),
                numpy_rng_state=rng.bit_generator.state, reverse_numpy_rng_state=reverse_rng.bit_generator.state,
                torch_rng_state=torch.get_rng_state(), generator_state=noise.get_state()), args.out/'LAST.pt')
            print(json.dumps(dict(event='validation', **row, selected_update=chosen_update)), flush=True)

        evaluate(0)
        for update in range(1, release['updates']+1):
            density.model.train()
            optimizer.zero_grad(set_to_none=True)
            # Inclusive pairs span the entire declared TRAIN panel. Their GT
            # never enters the separately sampled reverse-proposal condition.
            ids = rng.integers(0, len(train['x']), size=release['inclusive_NLL_batch'])
            prep = density.prepare(train, ids)
            inclusive = -density.log_prob(density.tensor(train['x'][ids]), prep).mean()
            lam = lambda_at(update, release['lambda_schedule'])
            phi_mean = entropy = reverse_value = prior_mean = None
            loss = inclusive
            if args.mode == 'physics':
                ii = physics_ids[reverse_rng.integers(0,len(physics_ids),size=4)]
                # No TRAIN labels are supplied to prepare/reparameterize.
                prep_reverse = density.prepare(train, ii)
                eps = torch.randn((4,66),generator=noise).to(density.device)
                x, logq = density.reparameterize(eps, prep_reverse)
                prior = .5*x.square().sum(1)
                book.context = dict(update=update, role='training_reverse', parents=train['parent'][ii].tolist(),lambda_value=lam)
                raw=x.detach().cpu().double().numpy()
                lq=logq.detach().cpu().double().numpy()
                prior_saved=prior.detach().cpu().double().numpy()
                queries=np.arange(oracle.query_count+1,oracle.query_count+5)
                batch_folder=book.jobdir/f'transport_batches/training/u{update:05d}'
                persist_transport_batch(batch_folder,'BEFORE',raw,lq,prior_saved,train['parent'][ii],update,lam,queries)
                try:
                    phi = ActualPhi.apply(x, original_batch(observations,train['parent'][ii]), oracle)
                except BaseException as exc:
                    persist_transport_batch(batch_folder,'AFTER',raw,lq,prior_saved,train['parent'][ii],update,lam,queries,error=repr(exc))
                    raise
                persist_transport_batch(batch_folder,'AFTER',raw,lq,prior_saved,train['parent'][ii],update,lam,queries,
                                        phi=phi.detach().cpu().double().numpy())
                reverse = (phi+prior+logq).mean()
                loss = inclusive+lam*reverse
                phi_mean, prior_mean = float(phi.detach().mean().cpu()), float(prior.detach().mean().cpu())
                entropy, reverse_value = float(-logq.detach().mean().cpu()), float(reverse.detach().cpu())
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError('nonfinite training objective')
            loss.backward()
            norm = torch.nn.utils.clip_grad_norm_(density.model.parameters(), release['grad_clip'], error_if_nonfinite=True)
            optimizer.step()
            completed = update
            train_metrics = dict(total=float(loss.detach().cpu()), inclusive_NLL=float(inclusive.detach().cpu()),
                lambda_value=lam if args.mode=='physics' else 0., Phi=phi_mean, prior_energy=prior_mean,
                entropy=entropy, reverse_ELBO=reverse_value, gradient_norm_before_clip=float(norm.detach().cpu()))
            if update in plan['validation_updates']:
                evaluate(update)
        if costs['RHS'] != plan['planned_RHS']:
            raise AssertionError('completed fit differs from frozen full physics roster')
    except BaseException as exc:
        status = 'BUDGET_STOP' if isinstance(exc, BudgetStop) else 'FAILED'
        error = repr(exc)
        raise
    finally:
        write(args.out/'RESULT.json', contract | dict(status=status, error=error,
            selected_update=chosen_update, selected_validation_metric=chosen_value if chosen_state is not None else None,
            completed_updates=completed, history=history, costs=costs, plan=plan, source_checkpoint=release['checkpoint'],
            wall_seconds=time.monotonic()-started, cpu_seconds=time.process_time()-cpu_started,
            posterior_qualification='NOT_ESTABLISHED', truth_at_deployment=False,
            checkpoint=str(args.out/'BEST.pt') if chosen_state is not None else None, receipt=str(book.jobdir/'CLOSED.json')))
        book.close(status,error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--job', required=True)
    parser.add_argument('--mode', choices=['physics','supervised'], required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--device', choices=['cpu','mps'], required=True)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
