#!/usr/bin/env python3
"""Standardize posterior files and collapse three ancestries to two."""

import configparser
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TWO_STATES = ["2,0", "1,1", "0,2"]
THREE_STATES = ["2,0,0", "1,1,0", "1,0,1", "0,2,0", "0,1,1", "0,0,2"]
DOSAGE = {s: tuple(map(int, s.split(","))) for s in THREE_STATES}


def read_config(path):
    c = configparser.ConfigParser()
    c.read(path)
    return c


def read_file_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip() and not x.startswith("#")]


def read_chromosomes(path):
    allowed = set()
    with open(path) as f:
        for line in f:
            p = line.strip().split()
            if not p or p[0].lower() in {"chrom", "chr", "chromosome"}:
                continue
            allowed.add(str(p[0]))
    if not allowed:
        raise ValueError(f"No chromosomes found in {path}")
    return allowed


def contribution_three(df):
    p = df[THREE_STATES].astype(float).fillna(0.0).to_numpy()
    dosage = np.array([DOSAGE[s] for s in THREE_STATES], dtype=float)
    site_total = p.sum(axis=1)
    expected = np.divide(
        p @ dosage,
        site_total[:, None],
        out=np.zeros((len(df), 3)),
        where=site_total[:, None] > 0,
    )
    return np.nanmean(expected / 2.0, axis=0)


def contribution_two(df):
    p = df[TWO_STATES].astype(float).fillna(0.0)
    total = p.sum(axis=1)
    major = (2.0 * p["2,0"] + p["1,1"]) / (2.0 * total.replace(0, np.nan))
    major = major.fillna(0.0).mean()
    return np.array([major, 1.0 - major])


def collapse_three(df, rank_order, third_state_max, apply_third_filter):
    p = df[THREE_STATES].astype(float).fillna(0.0)
    major, minor, third = rank_order
    out = np.zeros((len(df), 3), dtype=float)

    for state in THREE_STATES:
        a, b, c = DOSAGE[state]
        retained = (a if major == 0 else b if major == 1 else c,
                    a if minor == 0 else b if minor == 1 else c)
        if retained == (2, 0):
            col = 0
        elif retained == (1, 1):
            col = 1
        elif retained == (0, 2):
            col = 2
        else:
            continue
        out[:, col] += p[state].to_numpy()

    if apply_third_filter and third_state_max < 1.0:
        third_hom = [
            s for s in THREE_STATES
            if DOSAGE[s][third] == 2 and DOSAGE[s][major] == 0 and DOSAGE[s][minor] == 0
        ]
        if third_hom:
            mask = p[third_hom[0]].to_numpy() <= third_state_max
            out[~mask, :] = 0.0

    total = out.sum(axis=1)
    out = np.divide(out, total[:, None], out=np.zeros_like(out), where=total[:, None] > 0)
    return out


def process_file(path, output_path, method, allowed_chroms, third_state_max):
    df = pd.read_csv(path, sep=r"\s+", dtype={"chrom": str})
    if len(df.columns) < 3:
        raise ValueError(f"{path}: posterior file has too few columns")

    chrom_col = "chrom" if "chrom" in df.columns else df.columns[0]
    pos_col = "pos" if "pos" in df.columns else "position" if "position" in df.columns else df.columns[1]
    df[chrom_col] = df[chrom_col].astype(str)
    df = df[df[chrom_col].isin(allowed_chroms)].copy()
    if df.empty:
        raise ValueError(f"{path}: no sites remain after chromosome filtering")

    state_columns = set(df.columns[2:])
    if set(THREE_STATES).issubset(state_columns):
        ancestry_count = 3
        contribution = contribution_three(df)
        ranking = list(np.argsort(-contribution))
        fixed = collapse_three(
            df,
            ranking,
            third_state_max,
            apply_third_filter=(method == "LKsite"),
        )
        ranking_text = "THREE_TO_TWO"
    elif set(TWO_STATES).issubset(state_columns):
        ancestry_count = 2
        contribution = contribution_two(df)
        ranking = list(np.argsort(-contribution))
        fixed = df[TWO_STATES].astype(float).fillna(0.0).to_numpy()
        totals = fixed.sum(axis=1)
        fixed = np.divide(fixed, totals[:, None], out=np.zeros_like(fixed), where=totals[:, None] > 0)
        ranking_text = "TWO_ANCESTRY"
    else:
        raise ValueError(f"{path}: expected two- or three-ancestry genotype-state columns")

    sample = Path(path).name.split(".")[0]
    result = pd.DataFrame({
        "chrom": df[chrom_col].astype(str),
        "pos": pd.to_numeric(df[pos_col], errors="raise").astype(int),
        "2,0": fixed[:, 0],
        "1,1": fixed[:, 1],
        "0,2": fixed[:, 2],
    })
    result.to_csv(output_path, sep="\t", index=False)

    ranking_path = output_path.with_suffix(output_path.suffix + ".ranking.txt")
    with open(ranking_path, "w") as f:
        f.write(f"ANCESTRY_COUNT\t{ancestry_count}\n")
        f.write(f"STATUS\t{ranking_text}\n")
        for rank, idx in enumerate(ranking, 1):
            f.write(f"Rank{rank}_original\t{idx + 1}\t{contribution[idx]:.10f}\n")

    print(f"{Path(path).name}: ancestry={ancestry_count}, ranking={np.array(ranking) + 1}, contribution={contribution}")
    return ancestry_count, contribution, ranking


def main(config_path):
    c = read_config(config_path)
    root = Path(config_path).resolve().parent
    file_list = root / c["INPUT"]["FILE_LIST"]
    chrom_file = root / c["INPUT"]["CHROM_FILE"]
    outdir = root / c["OUTPUT"]["PREPARED_DIR"]
    outdir.mkdir(parents=True, exist_ok=True)
    metadata_dir = root / c["OUTPUT"]["BASELINE_DIR"]
    metadata_dir.mkdir(parents=True, exist_ok=True)

    method = c["GENERAL"]["METHOD"].strip()
    if method not in {"LKhard", "LKsite"}:
        raise ValueError("METHOD must be LKhard or LKsite")

    allowed = read_chromosomes(chrom_file)
    third_max = float(c["PARAMETER"]["LKSITE_THIRD_STATE_MAX"])

    summary = []
    for path_str in read_file_list(file_list):
        path = Path(path_str)
        if not path.is_absolute():
            path = (root / path).resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        out = outdir / path.name
        ancestry_count, contribution, ranking = process_file(path, out, method, allowed, third_max)
        if ancestry_count == 3:
            h = float(contribution[ranking[1]])
        else:
            h = float(min(contribution[0], contribution[1]))
        summary.append({"ID": path.name.split(".")[0], "h_minor_parental": h})

    pd.DataFrame(summary).to_csv(metadata_dir / "ancestry_h.txt", sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python 00.prepare_posterior.py config.ini")
    main(sys.argv[1])
