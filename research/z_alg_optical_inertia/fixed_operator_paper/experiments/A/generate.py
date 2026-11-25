"""Create the prospective observations and seal, after P0 and before any decision."""
from datetime import datetime, timezone
import secrets

import numpy as np
from scipy.linalg import solve, svdvals
from scipy.special import hankel1, jv

from common import BGS, PUBLIC, PRIVATE, RATIOS, REGION, RESULTS, SEAL, check_preregistration
from common import digest, forward, require_p0, save, write_json


def model():
    side = 5
    step = 0.7 / side
    k = 2 * np.pi
    a = step / np.sqrt(np.pi)
    y, x = np.indices((side, side))
    xy = (np.column_stack((x.ravel(), y.ravel())) - (side - 1) / 2) * step
    d = np.linalg.norm(xy[:, None] - xy[None, :], axis=2)
    d[d == 0] = 1
    coeff = 0.5j * np.pi * k * a * jv(1, k * a)
    green = coeff * hankel1(0, k * d)
    np.fill_diagonal(green, 0.5j * np.pi * k * a * hankel1(1, k * a) - 1)
    rx = -np.pi / 2 + np.arange(12) * np.pi / 12
    tx = np.arange(8) * np.pi / 4
    receive_xy = 2 * np.column_stack((np.cos(rx), np.sin(rx)))
    sensor = coeff * hankel1(0, k * np.linalg.norm(receive_xy[:, None] - xy[None, :], axis=-1))
    incident = np.exp(1j * k * xy @ np.vstack((np.cos(tx), np.sin(tx))))
    return dict(G=green, S=sensor, B=incident, V=np.eye(len(tx)), xy=xy,
                receiver_xy=receive_xy, receiver_angles=rx, source_angles=tx)


SHAPES = (
    (0, 1, 2, 3, 5),              # contiguous solid patch
    (0, 1, 2, 3, 4, 5, 6, 7),     # discrete annulus
    (0, 3, 4, 7),                 # two disconnected components
    (0, 1, 3, 5, 6, 7),           # graded open patch
)


def orient(support, rotation, flip):
    mask = np.zeros((5, 5), dtype=bool)
    mask.ravel()[REGION[list(support)]] = True
    mask = np.rot90(mask, rotation)
    if flip:
        mask = np.fliplr(mask)
    return np.flatnonzero(mask)


def ref_width(green, center, support):
    system = np.eye(len(center)) - center[:, None] * green
    field_map = solve(system.T, green.T, check_finite=False).T
    return 1 / svdvals(field_map[:, support], check_finite=False)[0]


def main():
    prereg = check_preregistration()
    p0_digest = require_p0()
    if SEAL.exists() or PRIVATE.exists() or PUBLIC.exists():
        raise FileExistsError("Prospective observations already exist")
    RESULTS.mkdir(exist_ok=True, mode=0o700)
    seed = secrets.token_bytes(16)
    rng = np.random.default_rng(int.from_bytes(seed))
    geom = model()
    green, sensor, incident = (geom[name] for name in ("G", "S", "B"))
    n = len(green)
    outside = np.setdiff1d(np.arange(n), REGION)
    centers = np.zeros((24, 2, n))
    widths = np.zeros((24, 2))
    bounds = np.zeros((24, len(BGS), len(RATIOS), 2, 2, n))
    truth = np.zeros((24, 2, len(BGS), n))
    data = np.zeros((24, 2, len(BGS), *sensor.shape[:1], incident.shape[1]), complex)
    clean = np.zeros_like(data)
    radii = np.zeros(data.shape[:3])
    branch = np.zeros((24, 2), dtype=int)
    pair_ids, target_ids, families = [], [], []
    for pair in range(24):
        amplitude = pair < 16
        shape = pair % 4
        rot, flip = rng.integers(4), rng.integers(2)
        other = shape if amplitude else (shape + 1 + (pair // 4) % 3) % 4
        supports = [orient(i, rot, flip) for i in (SHAPES[shape], SHAPES[other])]
        if amplitude:
            assert np.array_equal(*supports)
            peaks = [rng.uniform(.65, 2.2), rng.uniform(3.3, 7.5)]
        else:
            assert not np.array_equal(*supports)
            peaks = rng.uniform(1.8, 5.2, size=2)
        texture = rng.uniform(.60 if shape == 3 else .87, 1, n)
        for slot, support in enumerate(supports):
            pattern = texture[support]
            centers[pair, slot, support] = peaks[slot] * pattern / max(pattern)
            widths[pair, slot] = ref_width(green, centers[pair, slot], support)
        template_jitter = rng.uniform(-1, 1, size=n) * .2 * np.min(widths[pair]) / max(peaks)
        template_jitter = np.clip(template_jitter, -.10, .10)
        bg_texture = rng.uniform(.15, .85, size=len(outside))
        branch[pair] = rng.permutation(2)
        pair_ids.append(f"phase-a-pair-{pair+1:02d}")
        families.append("amplitude" if amplitude else "structure")
        target_ids.append([rng.bytes(8).hex(), rng.bytes(8).hex()])
        for bg, cap in enumerate(BGS):
            for w, ratio in enumerate(RATIOS):
                for slot, support in enumerate(supports):
                    center = centers[pair, slot, support]
                    width = widths[pair, slot] * ratio
                    box = bounds[pair, bg, w, slot]
                    box[0, support] = np.maximum(.015625, center - width)
                    box[1, support] = np.minimum(32., center + width)
                    box[1, outside] = cap
        for target, slot in enumerate(branch[pair]):
            noise = rng.standard_normal(data.shape[-2:]) + 1j * rng.standard_normal(data.shape[-2:])
            noise /= np.linalg.norm(noise)
            for bg, cap in enumerate(BGS):
                material = truth[pair, target, bg]
                material[:] = centers[pair, slot] * (1 + template_jitter)
                material[outside] = cap * bg_texture
                clean_block = forward(green, sensor, incident, material)[0]
                scale = np.linalg.norm(clean_block)
                clean[pair, target, bg] = clean_block
                data[pair, target, bg] = clean_block + .015 * scale * noise
                radii[pair, target, bg] = .015001 * scale
    save(PUBLIC, **geom, Y=data, radii=radii, centers=centers, ref_width=widths,
         boxes=bounds, pair_ids=np.array(pair_ids), target_ids=np.array(target_ids),
         families=np.array(families), bg=np.array(BGS), ratio=np.array(RATIOS))
    save(PRIVATE, material=truth, generating_branch=branch, clean=clean,
         seed=np.frombuffer(seed, dtype=np.uint8))
    PRIVATE.chmod(0o600)
    write_json(SEAL, dict(algorithm="BLAKE2b-256", files={p.name: digest(p) for p in (PUBLIC, PRIVATE)},
                          preregistration_digest=prereg["digest"], p0_digest=p0_digest,
                          freeze_utc=datetime.now(timezone.utc).isoformat(),
                          independent_pairs=24, target_materials=48, acquisition_blocks=192,
                          provenance="fresh confidential seed, no legacy objects or labels"))
    print(f"Frozen 24 pairs, {len(data.flat)} complex receiver scalars across all distinct blocks.",
          flush=True)


if __name__ == "__main__":
    main()
