"""Observation-only whole-child CPU budget hooks for inherited R8 inversion."""
import time


class CPUBudgetReached(BaseException):
    """Separate control flow: inherited except Exception must not label this failure."""


class WholeChildBudget:
    def __init__(self, seconds):
        self.seconds = float(seconds)
        if not 0 < self.seconds < float('inf'):
            raise ValueError('CPU budget must be positive and finite')
        self.checks = 0
        self.stop_cpu = None

    def check(self):
        self.checks += 1
        now = time.process_time()  # Absolute CPU since this child began, imports included.
        if now >= self.seconds:
            self.stop_cpu = now
            raise CPUBudgetReached()

    def observer_class(self, base):
        budget = self
        class BudgetObserver(base):
            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.eligible_best = None
                self.paid_states = []
                self.phase_costs = []
                self.incumbents = []
                self.last_paid_state = None
                budget.observer = self

            def record(self, *args, **kwargs):
                row = super().record(*args, **kwargs)
                self.last_paid_state = row
                completed = time.process_time()
                eligible = completed <= budget.seconds
                self.paid_states.append(dict(CPU_seconds=completed, phase=self.phase,
                    original_Phi=row['objective'], eligible_at_cap=eligible,
                    physical_model_cost=dict(self.model.cost)))
                if eligible and (self.eligible_best is None or row['objective'] < self.eligible_best['objective']):
                    self.eligible_best = row
                    self.incumbents.append((completed, row['m'].copy(), row['objective'], self.phase, dict(self.model.cost)))
                budget.check()
                return row

            def state(self, *args, **kwargs):
                budget.check()
                return super().state(*args, **kwargs)

            def physical_score(self, *args, **kwargs):
                budget.check()
                return super().physical_score(*args, **kwargs)

            def data_value_gradient(self, *args, **kwargs):
                budget.check()
                return super().data_value_gradient(*args, **kwargs)
        return BudgetObserver
