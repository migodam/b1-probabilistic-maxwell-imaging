"""Predesignated mature public-only controls with independent frozen budgets."""
from research.research_gn_nn_preservation_r4.src.anchors import generate_anchor

def run_baseline(operators,public,output,reference,kind):
    if kind not in ('fulltrf','weakreg1e-5','discrepancy'):raise ValueError(kind)
    return generate_anchor(operators,public,output,reference,1e-5 if kind=='weakreg1e-5' else 0.,dict(discrepancy_stop=kind=='discrepancy'))
