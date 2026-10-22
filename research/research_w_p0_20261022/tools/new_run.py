"""Prepare an explicitly separate reproduction; never resets prior budgets."""
from pathlib import Path
import sys,shutil
root=Path(__file__).resolve().parents[1]
dest=Path(sys.argv[1]).expanduser().resolve()
if dest.exists():raise RuntimeError('Destination must not exist; existing campaigns remain intact')
dest.mkdir(parents=True)
for name in ['src','tests','tools','source_package','docs']:
    shutil.copytree(root/name,dest/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
for name in ['w_config.json','PHYSICS_PROTOCOL.md','REPRODUCE.md','requirements.txt']:
    shutil.copy2(root/name,dest/name)
print('Independent reproduction directory prepared; zero physical calls:',dest)
