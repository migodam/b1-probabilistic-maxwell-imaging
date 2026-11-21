"""Read-only remote benchmark. No files or settings are written on the remote host."""
import base64
import io
import json
import os
import subprocess
import sys
from pathlib import Path
import numpy as np

here=Path(__file__).resolve().parent
root=here.parents[1]
hist=root.parent/'limited_aperture_routed_dr_20261120'
name='la180_austria_rot90_eps3p5'
ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8']
if os.getenv('SRAP_SSH_KEY'):
    ssh += ['-i',os.environ['SRAP_SSH_KEY']]
ssh.append(os.environ['SRAP_SSH_HOST'])
with np.load(hist/'.runs/public'/f'{name}.npz') as f:
    pack={k:f[k] for k in ('E','S','Y','centers','k','radius')}
with np.load(here.parent/'checkpoints'/f'{name}__shared.npz') as f:
    pack['material']=f['before_m']
stream=io.BytesIO()
np.savez_compressed(stream,**pack)
import threadpoolctl
portable=base64.b64encode(Path(threadpoolctl.__file__).read_bytes()).decode()
bootstrap="import base64,sys,types\ntp=types.ModuleType('threadpoolctl');tp.__file__='<embedded-threadpoolctl>';sys.modules['threadpoolctl']=tp;exec(compile(base64.b64decode("+repr(portable)+"),tp.__file__,'exec'),tp.__dict__)\n"
source=bootstrap+"PAYLOAD="+repr(base64.b64encode(stream.getvalue()).decode())+'\n'+(here/'fixed_work.py').read_text()
local=subprocess.run([sys.executable,'-B','-'],input=source,text=True,capture_output=True,cwd=root,check=True)
(here/'mac_fixed_work.json').write_text(local.stdout)
cmd=r'cd /d D:\AI\SRAP_CIE_20261121 && D:\python\python.exe -B -'
remote=subprocess.run(ssh+[cmd],input=source,text=True,capture_output=True)
if remote.returncode: raise RuntimeError(remote.stderr[-1500:])
(here/'xinan_fixed_work.json').write_text(remote.stdout)
