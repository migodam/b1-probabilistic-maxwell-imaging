"""Algebraic reciprocity diagnostic for the archived Manitoba session.

The published VNA/switch paths are direction dependent. Raw transpose projection
is an analysis sensitivity only, never an estimate of an antenna-plane error.
"""
import argparse
import csv
import json
from pathlib import Path
import statistics
from time import perf_counter, process_time
from zipfile import ZipFile

import numpy as np
from flint import acb, arb, ctx

HERE = Path(__file__).resolve().parent
DEFAULT_ARCHIVE = (
    HERE.parents[3] / ".artifacts/z_alg_validation/manitoba/Tomographic_Dataset_1.zip"
    if len(HERE.parents) > 3 else Path("Tomographic_Dataset_1.zip")
)
SOURCE_PAPER = "https://home.cc.umanitoba.ca/~lovetrij/Research/PDFs/Journals/JP50.pdf"
FREQUENCIES = ("3", "3.5", "4", "4.5", "5", "5.5", "6")
EMPTY = ("tomographic_incident.txt", "tomographic_incident_2.txt")
TARGETS = {
    "nylon_large": ("tomographic_nylon_4p0_inch_diam.txt","tomographic_nylon_4p0_inch_calibrated.txt"),
    "nylon_small": ("tomographic_nylon_1p5_inch_diam.txt","tomographic_nylon_1p5_inch_calibrated.txt"),
    "nylon_pair": ("tomographic_nylon_cysl_0p5cm.txt","tomographic_nylon_cyls_0p5cm_calibrated.txt"),
    "polyethylene_e": ("tomographic_e_phantom.txt","tomographic_e_phantom_calibrated.txt"),
}
COUNTS = (6, 8, 12, 16, 24)
METRICS = (
    "raw_complex_projection_relative_mid_scattering",
    "raw_complex_projection_relative_first_empty_scattering",
    "raw_complex_projection_relative_last_empty_scattering",
    "raw_phase_only_optimistic_floor_relative_mid_scattering",
    "supplied_Ez_complex_projection_relative",
    "minimum_common_empty_radius_relative_mid_scattering",
    "raw_loaded_complex_projection_relative_loaded_norm",
)


def load(archive, member):
    entries = {}
    header = []
    for line in archive.read(member).decode("latin-1").splitlines():
        if line.startswith("//"):
            header.append(line[2:].strip())
            continue
        if not line.strip():
            continue
        frequency, a, b, angle, elevation, real, imag = line.split()
        f = float(frequency)
        tx, rx = int(a), int(b)
        key = f, tx, rx
        if (key in entries or tx == rx or tx not in range(1,25) or rx not in range(1,25)
                or f not in map(float,FREQUENCIES) or float(angle) or float(elevation)):
            raise ValueError(f"Unexpected or duplicate antenna record in {member}: {key}")
        # Keep the original stored sign. A common conjugation would leave every
        # reported Frobenius distance and magnitude unchanged.
        entries[key] = acb(arb(real),arb(imag))
    if len(entries) != 7*24*23:
        raise ValueError(f"Expected all off-diagonal directions, got {len(entries)} in {member}")
    return entries, header


def square(z):
    return z.real*z.real + z.imag*z.imag


def matrix_pair_distance(data, direction, phase_only=False):
    # For one reciprocal pair a,b the closest equal value is their mean.
    # Count each pair once. Two squared directional residuals give 1/2 |a-b|².
    if phase_only:
        return (sum(((abs(data[a,b])-abs(data[b,a]))**2 for a,b in direction),arb(0))/2).sqrt()
    return (sum((square(data[a,b]-data[b,a]) for a,b in direction),arb(0))/2).sqrt()


def norm(data):
    return sum((square(x) for x in data.values()),arb(0)).sqrt()


def enclosure(value):
    if not value.is_finite():
        raise ArithmeticError("Non-finite interval")
    lower = np.nextafter(float(value.lower()), -np.inf)
    upper = np.nextafter(float(value.upper()), np.inf)
    return [float(lower),float(upper)]


def summarize(rows):
    result = {}
    for count in COUNTS:
        group = [r for r in rows if r["ports_each_side"]==count]
        components = {}
        for k in METRICS:
            lo,hi = zip(*(r[k] for r in group),strict=True)
            components[k] = dict(min_lower=min(lo),median_display=statistics.median((a+b)/2 for a,b in zip(lo,hi)),
                                 max_upper=max(hi),strictly_above_3pct=sum(v > .03 for v in lo))
        result[count] = dict(independent_objects=4,paired_frequency_rows=len(group),observed_directions=count*(count-1),
                             reciprocal_pairs=count*(count-1)//2,metrics=components)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive",type=Path,default=DEFAULT_ARCHIVE)
    args = parser.parse_args()
    start = perf_counter(),process_time()
    ctx.prec=128
    members = list(EMPTY)+[path for files in TARGETS.values() for path in files]
    with ZipFile(args.archive) as archive:
        entries={member:load(archive,member) for member in members}
    e0,e1=[entries[path][0] for path in EMPTY]
    rows=[]
    for q in COUNTS:
        subset=1+np.arange(q)*24//q
        pairs=[(int(a),int(b)) for a in subset for b in subset if a<b]
        directions=pairs+[(b,a) for a,b in pairs]
        for target,(tp,cp) in TARGETS.items():
            raw,calibrated=entries[tp][0],entries[cp][0]
            for freq in map(float,FREQUENCIES):
                loaded={(a,b):raw[freq,a,b] for a,b in directions}
                u={(a,b):calibrated[freq,a,b] for a,b in directions}
                delta_first={(a,b):loaded[a,b]-e0[freq,a,b] for a,b in directions}
                delta_last={(a,b):loaded[a,b]-e1[freq,a,b] for a,b in directions}
                delta_center={(a,b):(delta_first[a,b]+delta_last[a,b])/2 for a,b in directions}
                drift={(a,b):e1[freq,a,b]-e0[freq,a,b] for a,b in directions}
                amplitude=norm(delta_center)
                if not amplitude>0:
                    raise ArithmeticError("No direct scattering contrast in this slice")
                metrics={
                    METRICS[0]:matrix_pair_distance(delta_center,pairs)/amplitude,
                    METRICS[1]:matrix_pair_distance(delta_first,pairs)/norm(delta_first),
                    METRICS[2]:matrix_pair_distance(delta_last,pairs)/norm(delta_last),
                    METRICS[3]:matrix_pair_distance(delta_center,pairs,phase_only=True)/amplitude,
                    METRICS[4]:matrix_pair_distance(u,pairs)/norm(u),
                    METRICS[5]:norm(drift)/(2*amplitude),
                    METRICS[6]:matrix_pair_distance(loaded,pairs)/norm(loaded),
                }
                rows.append(dict(ports_each_side=q,target=target,frequency_GHz=freq,transmission_directions=len(directions),
                                 unordered_reciprocal_pairs=len(pairs),same_ports=[int(i) for i in subset],
                                 raw_scattered_mid_Fro=enclosure(amplitude),
                                 raw_loaded_Fro=enclosure(norm(loaded)),
                                 raw_two_empty_recipe_difference_Fro=enclosure(norm(drift)),
                                 raw_distance_to_literal_symmetric_matrix_Fro=enclosure(matrix_pair_distance(delta_center,pairs)),
                                 **{k:enclosure(v) for k,v in metrics.items()}))
    report=dict(source_paper=f"{SOURCE_PAPER}#page=4", acquisition="2010-07-28 Manitoba Dataset 1",
                raw_archive=args.archive.name, archive_members={key:entries[key][1] for key in entries},
                matrix_contract="Original VNA reference-plane transmissions, one direction-dependent switched path, off-diagonal only",
                supplied_fields_contract="Author's pairwise model-calibrated E_z, not common-plane scattering amplitudes",
                half_repeat="min common radius of two empty readings, or of two raw loaded-minus-empty recipes using the same loaded record",
                limitations="No repeated loaded target, no measured additive antenna impedance ΔZ, no antenna-plane reciprocity test, no 3% physical-noise claim",
                normalization="Each asymmetry uses the norm of the same off-diagonal selected data on both directions",
                computed_scope="Pure algebraic sensitivity in native file units, phases retained and separately tested as freely adjustable",
                arb_precision_bits=ctx.prec,summary=summarize(rows),
                elapsed_wall=perf_counter()-start[0],elapsed_CPU=process_time()-start[1])
    out=HERE/"results"
    (out/"manitoba_chamber_reciprocity.json").write_text(json.dumps(report,indent=2)+"\n")
    csv_rows=[]
    for row in rows:
        vals={}
        for k,v in row.items():
            if k=="same_ports":
                vals[k]=",".join(map(str,v))
            elif isinstance(v,list):
                vals[k+"_lower"],vals[k+"_upper"]=v
            else:
                vals[k]=v
        csv_rows.append(vals)
    with (out/"manitoba_chamber_reciprocity.csv").open("w") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    for q,count in report["summary"].items():
        metrics=count["metrics"]
        print(f"{q:2d} ports: {count['observed_directions']} observed directions, "
              f"{count['paired_frequency_rows']} target-frequency rows. "
              f"Median raw projection {100*metrics[METRICS[0]]['median_display']:.2f}%, "
              f"half empty-repeat radius {100*metrics[METRICS[5]]['median_display']:.2f}%.")


if __name__=="__main__":
    main()
