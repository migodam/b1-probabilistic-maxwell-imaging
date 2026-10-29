"""Freeze completed predictions before independent validation or truth access.

No digests are computed. Each manifest identifies immutable snapshot paths by
byte length and nanosecond modification time, and captures the method config.
"""
import argparse
import datetime
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pattern', required=True, action='append')
    p.add_argument('--name', required=True)
    a = p.parse_args()
    destination = ROOT / 'frozen_predictions' / a.name
    if destination.exists():
        raise FileExistsError(destination)
    sources = sorted({source for pattern in a.pattern for source in ROOT.glob(pattern)})
    if not sources:
        raise ValueError('No predictions matched')
    rows = []
    destination.mkdir(parents=True)
    for source in sources:
        result = json.loads(source.read_text())
        if 'candidates' not in result or 'accounting' not in result:
            raise ValueError('Not a prediction receipt: ' + str(source))
        target = destination / source.name
        shutil.copyfile(source, target)
        target.chmod(0o444)
        stat = target.stat()
        rows.append(dict(source=str(source.relative_to(ROOT)),
                         snapshot=str(target.relative_to(ROOT)),
                         bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                         method=result['method'], config=result['config'],
                         candidates=len(result['candidates']),
                         status=result['status']))
    manifest = dict(schema='R3_PREDICTION_FREEZE_V1',
                    frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    files=rows, truth_read=False, solver_feedback_after_freeze=False,
                    identity='Immutable snapshot plus byte length/mtime_ns; no SHA checks')
    path = ROOT / 'audit' / (a.name + '_FREEZE.json')
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(manifest, indent=2) + '\n')
    path.chmod(0o444)
    print(json.dumps(dict(manifest=str(path), predictions=len(rows), truth_read=False)))


if __name__ == '__main__':
    main()
