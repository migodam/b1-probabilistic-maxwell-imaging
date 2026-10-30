"""R4 uses the frozen R3 physical model without a second implementation."""
from research.research_z_alg_r3_digits.src.physics import Model, build_operators
from research.research_z_alg_r3_digits.src.multiscale import Regularizer, material_basis, resample_material

__all__ = ['Model', 'build_operators', 'Regularizer', 'material_basis', 'resample_material']
