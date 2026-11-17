from common import *
import platform,importlib.util,datetime
files=list((ROOT/'src').glob('*.py'))+[ROOT/'EXPERIMENT_PROTOCOL.md']
save(ROOT/'inputs/solver_freeze.json',dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files={str(f.relative_to(ROOT)):digest(f) for f in files},historical_files={str(f.relative_to(ROOT.parent)):digest(f) for f in (R2/'src').glob('*.py')},python=sys.version,platform=platform.platform(),homotopy=dict(status='NOT_RUN',phcpy_available=importlib.util.find_spec('phcpy') is not None,reason='No installed PHCpack or Julia executable; historical homegrown tracker is not a mature independent control. No new backend installed.')))
print('Solver and design frozen before generation.')
