import argparse,sys
from pathlib import Path
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S.parents[1]))
from research.research_gn_nn_preservation_r4.src.plots import make_plots
p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,default=S/'results/TEST_SCORES.json');p.add_argument('--runs',nargs='*');a=p.parse_args()
print(make_plots(a.scores,S,a.runs))
