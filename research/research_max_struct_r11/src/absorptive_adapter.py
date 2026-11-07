"""Process-local adapter over frozen BothBlocks current/PR/FFT runner, no source edits."""
from pathlib import Path
import json
from research.research_max_struct_r11.src.cie_feedback_current import run_cie_fft as frozen_run

def run_cie_fft(model,Y,m0,config):
    op=config['matrix_denominator'];label='absorptiveGreenBothBlocks' if op.rank else 'absorptiveScalarBothBlocks'
    meta=dict(current_transform=label,material_transform=label,denominator_shift='c=-1j*beta; outer prefactor real beta',auxiliary_residual_metric='Absorptive full moving denominator, original Maxwell numerator and physics',matrix_rank=op.rank,
      scalar_control_scope='K0 and same8-step matrix readout; real-CIE seed is only a proposal, not mature original CIE',novelty_scope='Fixed DEV mechanism audit; imaginary shift itself is not claimed novel')
    try:
        result=frozen_run(model,Y,m0,config)
        for row in result.get('history',[]):
            if 'iteration' in row:row.update(current_transform=label,material_transform=label)
        result['execution_receipt'].update(meta)
        return result
    finally:
        path=Path(config['receipt_path'])
        if path.exists():
            receipt=json.loads(path.read_text());receipt.update(meta);path.write_text(json.dumps(receipt,indent=2)+'\n')
