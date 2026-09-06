#!/usr/bin/env python3
"""Calculate observed or randomized LKhard/LKsite scores."""

import configparser
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

STATES = ["2,0", "1,1", "0,2"]
STATE_INDEX = {s: i for i, s in enumerate(STATES)}


def read_config(path):
    c = configparser.ConfigParser()
    c.read(path)
    return c


def sample_id(name):
    if ".rand.all." in name:
        return name.split(".rand.all.")[0]
    if ".rand." in name:
        return name.split(".rand.")[0]
    return name.split(".fixed.posterior.txt")[0].split("_ahmm2bed")[0].split(".posterior")[0]


def calculate_hwe(prepared_files):
    rows = []
    for f in prepared_files:
        df = pd.read_csv(f, sep="\t")
        s20, s11, s02 = [df[s].astype(float).sum() for s in STATES]
        total = s20 + s11 + s02
        if total <= 0:
            continue
        a = (s20 + 0.5 * s11) / total
        b = 1.0 - a
        rows.append({"ID": sample_id(f.name), "P20": a * a, "P11": 2 * a * b, "P02": b * b})
    return pd.DataFrame(rows)


def load_hwe(path):
    df = pd.read_csv(path, sep="\t")
    return df.set_index("ID").to_dict("index")


def calculate_lkhard_segments(files, hwe):
    lookup = load_hwe(hwe) if isinstance(hwe, (str, Path)) else hwe.set_index("ID").to_dict("index")
    all_breakpoints = set()
    samples = []

    for f in files:
        sid = sample_id(f.name)
        if sid not in lookup:
            raise ValueError(f"{sid}: HWE probabilities not found")
        segments = []
        with open(f) as g:
            next(g, None)
            for line in g:
                p = line.rstrip().split("\t")
                if len(p) < 4:
                    continue
                try:
                    chrom, start, end, state = p[0], int(p[1]), int(p[2]), p[3]
                except ValueError:
                    continue
                if end > start:
                    segments.append((chrom, start, end, state))
                    all_breakpoints.add((chrom, start))
                    all_breakpoints.add((chrom, end))
        samples.append((sid, segments, lookup[sid]))

    if not samples:
        raise RuntimeError("No valid AHMM files found")

    by_chrom = {}
    for chrom, pos in all_breakpoints:
        by_chrom.setdefault(chrom, set()).add(pos)

    rows = []
    for chrom in sorted(by_chrom, key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x))):
        bp = sorted(by_chrom[chrom])
        for start, end in zip(bp[:-1], bp[1:]):
            if end <= start:
                continue
            likelihood = 1.0
            informative = False
            for sid, segments, probs in samples:
                state = None
                for ch, s, e, st in segments:
                    if ch != chrom:
                        continue
                    if s <= start and e >= end:
                        state = st
                        break
                if state in STATE_INDEX:
                    informative = True
                    likelihood *= [probs["P20"], probs["P11"], probs["P02"]][STATE_INDEX[state]]
            if informative:
                score = -np.log10(max(likelihood, np.finfo(float).tiny))
                rows.append([chrom, start, end, likelihood, score])

    return pd.DataFrame(rows, columns=["chrom", "pos", "end", "likelihood", "score"])


def calculate_lksite(files, h_map):
    combined = []
    for f in files:
        sid = sample_id(f.name)
        if sid not in h_map:
            raise ValueError(f"{sid}: h not found")
        h = h_map[sid]
        df = pd.read_csv(f, sep="\t", dtype={"chrom": str})
        one_minus_h = 1.0 - h
        df["individual_lk"] = (
            df["0,2"] * h ** 2
            + df["2,0"] * one_minus_h ** 2
            + df["1,1"] * h * one_minus_h
        )
        combined.append(df[["chrom", "pos", "individual_lk"]])
    if not combined:
        raise RuntimeError("No valid posterior files found")
    x = pd.concat(combined, ignore_index=True)
    result = x.groupby(["chrom", "pos"], as_index=False)["individual_lk"].prod()
    result = result.rename(columns={"individual_lk": "likelihood"})
    result["score"] = -np.log10(result["likelihood"].clip(lower=np.finfo(float).tiny))
    return result


def main(config_path, mode):
    c = read_config(config_path)
    root = Path(config_path).resolve().parent
    method = c["GENERAL"]["METHOD"].strip()
    prepared = root / c["OUTPUT"]["PREPARED_DIR"]
    hmm = root / c["OUTPUT"]["HMM_DIR"]
    random_dir = root / c["OUTPUT"]["RANDOM_DIR"]
    obsdir = root / c["OUTPUT"]["OBS_DIR"]
    base = root / c["OUTPUT"]["BASELINE_DIR"]
    obsdir.mkdir(parents=True, exist_ok=True)
    base.mkdir(parents=True, exist_ok=True)

    if method == "LKhard":
        prepared_files = sorted(prepared.glob("*.posterior"))
        hwe_file = base / "hwe_proportions_summary.txt"
        if mode == "observed":
            hwe = calculate_hwe(prepared_files)
            hwe.to_csv(hwe_file, sep="\t", index=False)
            files = sorted(hmm.glob("*_ahmm2bed"))
            result = calculate_lkhard_segments(files, hwe)
            out = obsdir / "observed_lkhard.txt"
        else:
            if not hwe_file.exists():
                raise FileNotFoundError(hwe_file)
            files = sorted(random_dir.glob("*.rand.*.fixed.bed"))
            result = calculate_lkhard_segments(files, hwe_file)
            out = base / "randomized_lkhard.txt"
    elif method == "LKsite":
        hfile = base / "ancestry_h.txt"
        if not hfile.exists():
            raise FileNotFoundError(hfile)
        hdf = pd.read_csv(hfile, sep="\t")
        h_map = dict(zip(hdf["ID"].astype(str), hdf["h_minor_parental"].astype(float)))
        if mode == "observed":
            files = sorted(prepared.glob("*.posterior"))
            result = calculate_lksite(files, h_map)
            out = obsdir / "observed_lksite.txt"
        else:
            files = sorted(random_dir.glob("*.rand.all.*.txt"))
            by_perm = {}
            for f in files:
                m = re.search(r"\.rand\.all\.(\d+)\.txt$", f.name)
                if m:
                    by_perm.setdefault(int(m.group(1)), []).append(f)
            outputs = []
            for k, group in sorted(by_perm.items()):
                r = calculate_lksite(group, h_map)
                r["permutation"] = k
                outputs.append(r)
            if not outputs:
                raise RuntimeError("No randomized LKsite files found")
            result = pd.concat(outputs, ignore_index=True)
            out = base / "randomized_lksite.txt"
    else:
        raise ValueError("METHOD must be LKhard or LKsite")

    result.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python 01.calculate_likelihood.py config.ini observed|randomized")
    main(sys.argv[1], sys.argv[2])
