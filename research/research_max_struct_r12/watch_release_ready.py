"""Finite readiness updates for remote03/04 while the final download continues."""
from pathlib import Path
import json,time
import recover_release_staging as recovery
for number in [3,4]:
    name=f'max-struct-r12-remote-{number:02d}.zip.zst'
    while not(recovery.OUT/name).exists():time.sleep(5)
    recovery.verify_available_remote()
    print(json.dumps(dict(accepted=name,scope='Verified available logical archive, not a publication action')),flush=True)
