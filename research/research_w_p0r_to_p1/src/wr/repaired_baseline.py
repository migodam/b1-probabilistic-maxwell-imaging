"""Frozen post-interface-repair L-BFGS-B replay, same BP/no GT initializer."""
from .online_inference import run_one
if __name__=='__main__':run_one(2026102201,'wide','LBFGSB_REPAIRED',iterations=18,max_calls_override=20)
