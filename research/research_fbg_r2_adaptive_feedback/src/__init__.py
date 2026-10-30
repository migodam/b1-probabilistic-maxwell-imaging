"""Adaptive fixed-feedback scalar-TM primitives; CUDA lane adapter lives elsewhere."""

from .feedback_fit import (FeedbackOperator, candidate_feedback, cap_feedback_at_anchors,
                           feedback_jvp, feedback_value_gradient, fit_feedback, fit_tangent_ridge,
                           source_normalized_current_tangents, spatial_2x2_indices)

__all__ = ["FeedbackOperator", "candidate_feedback", "cap_feedback_at_anchors",
           "feedback_jvp", "feedback_value_gradient", "fit_feedback", "fit_tangent_ridge",
           "source_normalized_current_tangents", "spatial_2x2_indices"]
