"""One stage-event control; no new metric, threshold, optimizer or observation."""
from dualmetric_core import StatefulDualMetric


class StageScheduledDualMetric(StatefulDualMetric):
    def _enter_stage(self,index):
        schedule=self.config['stage_metric_pairs']
        if len(schedule)!=len(self.stages):
            raise ValueError('Declare one metric pair per original stage')
        super()._enter_stage(index)
        # The inherited stage entry already resets PR and transports the same
        # coefficients. Only the next stage's two residual families change.
        self.switch_pair(str(schedule[index]),reset_memory=False)
