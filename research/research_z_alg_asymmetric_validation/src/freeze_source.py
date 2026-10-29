from pathlib import Path
import hashlib,json,platform,subprocess,time
from accounting import P,write

def main():
    files=[P/'EXPERIMENT_CONTRACT.json',*sorted((P/'src').glob('*.py'))]
    write(P/'SOURCE_MANIFEST.json',dict(status='FROZEN_BEFORE_NEW_PHYSICS',utc=time.time(),base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=P).decode().strip(),
        source_files=[dict(path=f.relative_to(P).as_posix(),bytes=f.stat().st_size,blake2b=hashlib.blake2b(f.read_bytes()).hexdigest()) for f in files],
        original_inputs='inputs/INPUT_IDENTITY.json CRC plus BLAKE2b; supplied SHA metadata retained without new SHA calculation',
        software=dict(python=platform.python_version(),platform=platform.platform()),
        amendments='New report/analysis helpers may be added; any core correction must preserve original source and explicit impact receipt before further physics. Frozen statistical/geometry criteria cannot change after outcomes.'))
if __name__=='__main__':main()
