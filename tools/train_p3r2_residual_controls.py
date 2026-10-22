#!/usr/bin/env python3
"""Separate legal-decoder controls; no change to original residual trainer defaults."""
import sys
from train_p3r2_residual import main
from p3r2.residual_risk import DEFAULT_CONTROL_MODELS


if __name__=='__main__':
    arguments=sys.argv[1:]
    if not any(arg=='--max-seconds' or arg.startswith('--max-seconds=') for arg in arguments):
        arguments+=['--max-seconds','360']
    main(argv=arguments,default_models=DEFAULT_CONTROL_MODELS)
