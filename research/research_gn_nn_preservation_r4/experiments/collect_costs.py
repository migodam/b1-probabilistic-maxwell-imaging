import sys
from pathlib import Path
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S.parents[1]))
from research.research_gn_nn_preservation_r4.src.costs import collect
print(collect(S))
