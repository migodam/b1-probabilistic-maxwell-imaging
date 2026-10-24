"""Prelocked revision1 exact bin labels. Parent launches all training jobs.

Reuses the frozen first-competition trainer mechanics and feature pipeline in a
separate process. Original risk.py/train_risk.py are never edited. The two
in-process class/label bindings select this revision's standalone mathematical
rule; checkpoint weights remain compatible with RiskActionPredictor.
"""
import argparse,json,sys
from pathlib import Path
import torch
from .risk_bin import BinRisk,bin_phase_labels,quantile_bins,MOMENT_RULE
from . import train_risk as original


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if '--help' in argv or '-h' in argv:
        print(__doc__);original.main(argv);return
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--kind',required=True,choices=('phase8','phase16'));parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--anchor-train','--train-context',required=True,type=Path)
    parser.add_argument('--anchor-validation','--validation-context',required=True,type=Path)
    args,_=parser.parse_known_args(argv)
    # These substitutions exist only in this dedicated process. Feature,
    # optimizer, split, selection, and resource reporting stay identical.
    original.LocalRisk=BinRisk;original.phase_labels=bin_phase_labels
    original.main(argv)
    checkpoint=args.out/'BEST.pt';p=torch.load(checkpoint,map_location='cpu',weights_only=False)
    k=int(args.kind[5:]);edges,centers=quantile_bins(k)
    metadata=dict(revision='PRELOCKED_RISK_REVISION1',moment_rule=MOMENT_RULE,
        exact_label='delta_component*clip(K*original_saturation(s)-j,0,1), j0..K-1; integration weights1/K',
        label_identity='sum_j label_j/K = delta*original_saturation(s) pointwise before float32 roundoff',
        label_semantics='h-prime weighted average joint event moment within equal h-mass quantile bin; not point-threshold event probability',
        bin_edges=edges.tolist(),bin_centers=centers.tolist(),bin_weights=[1/k]*k,
        correlated_property='same joint-parent delta and field; never E(delta)*E(event)',
        material_initialization_warning='uniform positive masses at zero logits: near half maximum legal increment; differs from anchor-initialized direct residual CNN',
        scientific_scope='original P1 decoder unchanged; removes finite point-quadrature label bias, not regression error or posterior qualification')
    p.update(metadata)
    # Atomic replacement only inside this new run's output directory.
    temporary=args.out/'BEST.bin_metadata.tmp';torch.save(p,temporary);temporary.replace(checkpoint)
    result_path=args.out/'RESULT.json';result=json.loads(result_path.read_text());result.update(metadata)
    result_path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (args.out/'BIN_RULE.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':'BIN_REVISION1_COMPLETE','out':str(args.out),'moment_rule':MOMENT_RULE,'parameter_count':result['parameter_count']}))

if __name__=='__main__':main()
