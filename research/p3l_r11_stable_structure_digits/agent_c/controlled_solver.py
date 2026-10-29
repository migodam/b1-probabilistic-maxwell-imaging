"""Inversion-compatible fixed total-field solver with a global loss gate.

Shared Galerkin mechanics are imported from Sol-B.  The global and local
stability decisions are separate, and every standalone instance pays setup.
"""
from __future__ import annotations
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'agent_b'))
sys.path.insert(0, str(HERE))
from fixed_total_field import TotalFieldSolver, CostLedger
from dissipative_bound import GlobalDissipativeBound


class _ForwardLedger(CostLedger):
    """Preserve detailed own counts and forward paid actions to a caller ledger."""
    def __init__(self, forward=None, path=None, rhs_cap=20000):
        self.forward = forward
        super().__init__(path=path, rhs_cap=rhs_cap)

    def event(self, event, **kw):
        super().event(event, **kw)
        if self.forward is not None:
            self.forward(dict(event=event, **kw))


class GlobalControlledSolver(TotalFieldSolver):
    name = 'GLOBAL_TOTAL_FIELD_CONTROLLED'

    def __init__(self, C, O, E, Fh, sigma=.00119211798514519,
                 ledger=None, *, a_anchor=None, ranks=(32,64,96,128),
                 basis_kind='green', uniform_port=False, rhs_cap=20000,
                 ledger_path=None, bound_mode='global'):
        if ledger is not None and hasattr(ledger, 'reserve') and hasattr(ledger, 'event'):
            own = ledger
        else:
            own = _ForwardLedger(ledger, path=ledger_path, rhs_cap=rhs_cap)
        gate = GlobalDissipativeBound(C, O, sigma, mode=bound_mode,
                                     uniform_port=uniform_port)
        self.bound_object = gate
        super().__init__(C, O, E, Fh, a_anchor=a_anchor, ranks=ranks,
                         basis_kind=basis_kind, ledger=own,
                         bound_callback=gate, sigma=sigma)

    def summary(self):
        return dict(name=self.name, basis=self.setup,
                    model_setup_seconds=self.model_setup_seconds,
                    bound=self.bound_object.summary(), ledger=self.ledger.summary(),
                    exact_derivatives=True,
                    reduced_derivative_status='NOT_CONTROLLED_NOT_USED')


ControlledGlobalSolver = GlobalControlledSolver
