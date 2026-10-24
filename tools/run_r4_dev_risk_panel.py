"""Execute authorized common exposedDEV panel only; zero training/Maxwell."""
from pathlib import Path
import subprocess,sys,json
ROOT=Path(__file__).resolve().parents[1];S=ROOT/'research/p3l_r4_conditional_transport';P=S/'results/development_risk_strong_panel_v1'
def run(*args):subprocess.run([sys.executable,str(ROOT/'tools/r4_risk_panel.py'),*map(str,args)],check=True,cwd=ROOT)
if __name__=='__main__':
 spec=json.loads((P/'SPEC.json').read_text())
 if spec['parents']!=[f'R4DEV_{i:05d}' for i in range(128)] or 'LOCKED' not in spec['cohort']:raise ValueError('fixed exposedDEV scope')
 run('freeze','--spec',P/'SPEC.json','--out',P/'FROZEN.json')
 run('features','--freeze',P/'FROZEN.json','--observations',S/'inputs/online','--out',P/'features')
 run('infer','--freeze',P/'FROZEN.json','--contexts',S/'results/contexts/dev','--native',S/'results/contexts/details','--online',P/'online','--out',P/'features')
 run('score','--freeze',P/'FROZEN.json','--online',P/'online','--truth',S/'inputs/offline_only','--out',P/'scores')
