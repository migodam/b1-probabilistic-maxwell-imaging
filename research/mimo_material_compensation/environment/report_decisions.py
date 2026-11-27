"""Independent data-only reading of sealed decisions, including the forbidden score gap."""
from pathlib import Path
import json
import numpy as np
from passive_replay import fro, rd, lower, up

HERE = Path(__file__).resolve().parent


def decide(side, contract):
    path = HERE / f"low_energy_{side}"
    proposal = json.loads((path / "proposal.json").read_text())
    previous = json.loads((path / "decision.json").read_text())
    direction = np.load(path / "selector.npz")["direction"]
    assert fro(direction) <= 1
    out = {}
    with np.load(path / "readings.npz") as frames:
        for specimen in frames:
            frame = frames[specimen]
            frame = np.r_[frame.real.ravel(), frame.imag.ravel()]
            value = float(direction@frame)
            err = up(rd(len(frame))*fro(frame)*fro(direction))
            lo = lower(value-err)
            hi = float(np.nextafter(value+err, np.inf))
            possible = []
            if lo <= proposal["low_read_score_upper"]:
                possible.append("low")
            if hi >= proposal["high_read_score_lower"]:
                possible.append("high")
            answer = possible[0] if len(possible)==1 else None
            if previous["results"][specimen]["label"]!=answer:
                raise ValueError("Gap decision disagrees with saved physical decision")
            out[specimen] = dict(possible=possible, score_with_error=[lo,hi],
                                roi_interval=contract["roi_outer_"+answer] if answer else None)
    return out


if __name__=="__main__":
    limits = json.loads((HERE / "decision_limits.json").read_text())["records"]
    for row in limits:
        if row.get("operating_energy") is None:
            continue
        result = decide(row["side"], row)
        path=HERE/f'low_energy_{row["side"]}/outward_decision.json'
        path.write_text(json.dumps(result,indent=2)+"\n")
