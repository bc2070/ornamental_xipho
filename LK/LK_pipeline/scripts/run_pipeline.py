#!/usr/bin/env python3
"""Run the unified LKhard/LKsite pipeline."""

import configparser
import subprocess
import sys
from pathlib import Path


def run(cmd, cwd):
    print("\n$ " + " ".join(map(str, cmd)))
    subprocess.run(cmd, cwd=cwd, check=True)


def main(config_path):
    config = Path(config_path).resolve()
    root = config.parent
    c = configparser.ConfigParser(); c.read(config)
    method = c["GENERAL"]["METHOD"].strip()
    python = sys.executable
    s = root / "scripts"

    steps = [
        ([python, str(s / "00.prepare_posterior.py"), str(config)], "Posterior preprocessing"),
        ([python, str(s / "02.posterior_to_ahmm.py"), str(config)], "Posterior to AHMM blocks"),
        ([python, str(s / "03.merge_ancestry_info.py"), str(config)], "Merge ancestry information"),
        ([python, str(s / "01.calculate_likelihood.py"), str(config), "observed"], "Observed likelihood"),
        ([python, str(s / "04.randomization.py"), str(config)], "Randomization"),
        ([python, str(s / "01.calculate_likelihood.py"), str(config), "randomized"], "Randomized likelihood"),
        ([python, str(s / "05.baseline.py"), str(config)], "Baseline threshold"),
    ]

    print("=" * 70)
    print(f"Unified LK pipeline: {method}")
    print(f"Root: {root}")
    print("=" * 70)
    for cmd, label in steps:
        print(f"\n[STEP] {label}")
        run(cmd, root)

    plot_prefix = root / c["OUTPUT"]["FIGURE_DIR"] / f"Manhattan_{method}"
    (root / c["OUTPUT"]["FIGURE_DIR"]).mkdir(parents=True, exist_ok=True)
    run(["Rscript", str(s / "06.plot_manhattan.R"), str(config), method, str(plot_prefix)], root)

    print("\n" + "=" * 70)
    print(f"Pipeline completed successfully: {method}")
    print("=" * 70)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python scripts/run_pipeline.py config.ini")
    main(sys.argv[1])
