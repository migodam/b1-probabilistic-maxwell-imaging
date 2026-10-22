"""Read-only environment metadata; no credentials, identifiers, or physical work."""
from pathlib import Path
import importlib.metadata as m
import platform, os, json, subprocess, sys
root=Path(__file__).resolve().parents[1]
packages={p:m.version(p) for p in ('numpy','scipy','matplotlib','scikit-fem','triangle')}
def sysctl(name):
    try:return subprocess.check_output(['sysctl','-n',name],text=True).strip()
    except (OSError,subprocess.CalledProcessError):return None
row=dict(python=sys.version,platform=platform.platform(),machine=platform.machine(),
         cpu_model=sysctl('machdep.cpu.brand_string'),logical_cpus=os.cpu_count(),
         physical_memory_bytes=sysctl('hw.memsize'),packages=packages,backend='CPU only',
         blas_threads=1,no_GPU_used=True,source_commit='db7146cae2d74ba231577d42f225dbfc9a389268',
         dependency_install_scope='scikit-fem/Triangle in project .deps; not redistributed',
         packing='source-major complex (4,12), real whitening sigma=0.01; norm is not divided by channel count')
(root/'docs/ENVIRONMENT.json').write_text(json.dumps(row,indent=2))
print(json.dumps(row))
