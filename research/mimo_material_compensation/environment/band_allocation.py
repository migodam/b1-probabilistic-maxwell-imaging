"""Necessary new-band energy from the same material witness across both frequencies."""
from fractions import Fraction
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent

def outward(q, upper):
    x=float(q)
    if (upper and Fraction(x)<q) or (not upper and Fraction(x)>q):
        x=math.nextafter(x, math.inf if upper else -math.inf)
    return x

def evaluate(side,gap,stem):
    u=[]
    filenames=[]
    for freq in (400,460):
        name=f"kyp_verified_responses_{side}_{gap}_{stem}_{freq}.json"
        rec=json.loads((HERE/name).read_text())
        if not rec['feasible']:
            raise ArithmeticError(f"Unverified norm in {name}")
        filenames.append(name)
        u.append(Fraction(rec['upper_over_2r']))
    # Total energy is exactly one; old observations must already admit both materials.
    a,b=(v*v for v in u)
    peak=max(a,b)
    row=dict(side=side,core_difference=gap,material_pair=stem,
             verified_upper_each_band_over_2r=[float(x) for x in u],sources=filenames,
             no_policy_works_for_total_E_at_most=outward(1/peak,False))
    if max(a,b)<=1:
        row.update(unit_energy_result='all allocations fail for a common history')
    elif a<1<b:
        row.update(unit_energy_result='a minimum fraction of the unit energy must be placed at 460 MHz',
                   necessary_E460_strictly_greater_than=outward((1-a)/(b-a),False))
    else:
        row.update(unit_energy_result='this pair does not impose a positive universal lower allocation to 460 MHz')
    return row

if __name__=="__main__":
    designs=[(32,2.0,'phase_32_two_shells'), (96,2.0,'wave_96_2.0_1'),(96,1.0,'phase_96_two_shells'),
             (32,1.5,'joint_wave_32_1.5_s8'),(96,1.0,'joint_wave_96_1.0_s4')]
    out=dict(energy_contract='total available incident energy, one joint absolute Frobenius noise ball; same nondispersive material in both bands; 16 active ports and one scalar passive feedback per shot',
             records=[evaluate(*design) for design in designs])
    (HERE/'band_allocation.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in out['records']:
        print(row['material_pair'],row.get('necessary_E460_strictly_greater_than'),row['no_policy_works_for_total_E_at_most'])
