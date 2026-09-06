#!/usr/bin/env python3
"""Calculate the null threshold and classify observed scores."""

import configparser
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main(config_path):
    c = configparser.ConfigParser(); c.read(config_path)
    root = Path(config_path).resolve().parent
    method = c["GENERAL"]["METHOD"].strip()
    base = root / c["OUTPUT"]["BASELINE_DIR"]
    obsdir = root / c["OUTPUT"]["OBS_DIR"]
    base.mkdir(parents=True, exist_ok=True)
    q = float(c["PARAMETER"]["THRESHOLD_PERCENTILE"])

    rand = pd.read_csv(base / f"randomized_{method.lower()}.txt", sep="\t")
    obs = pd.read_csv(obsdir / f"observed_{method.lower()}.txt", sep="\t")
    if rand.empty or "score" not in rand.columns:
        raise RuntimeError("Randomized likelihood file contains no score values")
    threshold = float(np.percentile(rand["score"].dropna(), q))

    pd.DataFrame({"percentile": [q], "threshold_score": [threshold]}).to_csv(
        base / "baseline_threshold.txt", sep="\t", index=False
    )
    obs["above_threshold"] = obs["score"] >= threshold
    obs.to_csv(base / "observed_vs_baseline.txt", sep="\t", index=False)
    print(f"{method} baseline ({q}%) = {threshold}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python 05.baseline.py config.ini")
    main(sys.argv[1])
