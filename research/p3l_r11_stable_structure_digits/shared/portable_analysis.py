"""Lossless evidence-only analysis entry. Does not launch any Maxwell solve."""
import sys,types,importlib.util,runpy,argparse
from pathlib import Path
HERE=Path(__file__).resolve().parent;P=HERE.parent
spec=importlib.util.spec_from_file_location('p1_material.prior',P/'inputs/PRIOR_SOURCE.py')
pkg=types.ModuleType('p1_material');pkg.__path__=[];sys.modules['p1_material']=pkg
pr=importlib.util.module_from_spec(spec);sys.modules['p1_material.prior']=pr;spec.loader.exec_module(pr)
sys.path.insert(0,str(HERE));import core
import numpy as np,json

def load_local():
 z=dict(np.load(P/'inputs/FIXED_HOST_SETUP.npz'));z['E']=z['ER'];return z

def local_prior():return pr.MaterialPrior(json.loads((P/'inputs/PRIOR_CONFIG.json').read_text()),np.load(P/'inputs/PRIOR_POINTS.npz')['points'])
core.load_setup=load_local;core.prior=local_prior
ap=argparse.ArgumentParser();ap.add_argument('analysis',choices=['score','cost']);a=ap.parse_args()
runpy.run_path(str(HERE/('score_and_plot.py' if a.analysis=='score' else 'finalize_cost_and_diagnostics.py')),run_name='__main__')
