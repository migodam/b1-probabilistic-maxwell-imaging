"""Read-only saved-bundle checker; does not import/instantiate a solver."""
import argparse
import json

from research_t.handoff import validate_bundle


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bundle")
    args = parser.parse_args()
    result = validate_bundle(args.bundle)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    raise SystemExit(0 if result["valid_bundle_semantics"] else 1)
