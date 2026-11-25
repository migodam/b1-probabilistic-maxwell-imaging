"""Set new-shot energy from old-data class bounds; all controls get the same exposure."""
from pathlib import Path
import json, shutil, sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'beams'))
from beam_common import save,atomic_write,digest
from datetime import datetime,timezone

def main(side):
    src=HERE/f'confirm_{side}'; dest=HERE/f'power_{side}'; dest.mkdir(exist_ok=True)
    if (dest/'future_readings.npz').exists():raise ValueError('use a clean output directory')
    old_design=json.loads((src/'proposal.json').read_text())
    array=np.load(src/'proposal_arrays.npz')
    control_src=HERE.parent/'beams'/f'confirm_{side}'
    ref=np.load(control_src/'proposal_arrays.npz')
    for file in ('registration.json','old_readings.npz','public_classes.npz'):shutil.copyfile(src/file,dest/file)
    output={};rows={}
    for idx,row in old_design['rows'].items():
        d=row['cases']['indefinite'].get('numerical_distance_over_r',0.)
        e0=row['fluence']; cap=1.
        e=min(cap,e0*(2.15/d)**2) if d>0 else cap
        prefix=idx+'__';V=array[prefix+'actions'];refV=ref[prefix+'actions']
        policies=row['policies']+['uniform_angle','random_angle','fisher_angle','unconditioned_source_eigenbeam']
        actions=np.column_stack((V,refV[:,[0,1,2,4]]))*np.sqrt(e/e0)
        output[prefix+'actions']=actions
        for tag in ('positive','indefinite'):output[prefix+tag+'_forms']=array[prefix+tag+'_forms']
        rows[idx]=dict(row,policies=policies,fluence=e,initial_fluence=e0,
                       threshold_candidate=e0*(2/d)**2 if d>0 else None,
                       unit_energy_cap=cap,design_target_over_r=2.15)
    save(dest/'proposal_arrays.npz',**output)
    atomic_write(dest/'proposal.json',dict(rows=rows,physical_source=f'confirm_{side}',design='same selected class quadratic, target 2.15r with maximum energy of one first shot; same energy and receiver error for all arms; no future readings input',
               design_utc=datetime.now(timezone.utc).isoformat(),
               design_binary_digest=digest(np.frombuffer((dest/'proposal_arrays.npz').read_bytes(),dtype=np.uint8)),
               old_file_digest=digest(np.frombuffer((dest/'old_readings.npz').read_bytes(),dtype=np.uint8))))
if __name__=='__main__':main(int(sys.argv[1]))
