#!/usr/bin/env python3
"""Convert standardized posterior probabilities into ancestry blocks."""

import configparser
import sys
from pathlib import Path

import pandas as pd

STATES = ["2,0", "1,1", "0,2"]


def natural_key(x):
    s = str(x)
    return (0, int(s)) if s.isdigit() else (1, s)


def load_chromosomes(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.strip().split()
            if not p or p[0].lower() in {"chrom", "chr", "chromosome"}:
                continue
            if len(p) >= 3:
                d[str(p[0])] = int(float(p[2]))
            elif len(p) >= 2:
                d[str(p[0])] = int(float(p[1]))
    if not d:
        raise ValueError(f"No chromosome lengths found in {path}")
    return d


def convert(df, lengths, threshold, max_gap):
    rows = []
    for chrom, sub in df.groupby("chrom", sort=False):
        chrom = str(chrom)
        if chrom not in lengths:
            continue
        sub = sub.sort_values("pos").reset_index(drop=True)
        chrom_end = lengths[chrom]
        start = int(sub.loc[0, "pos"])
        last_pos = start
        last_state = None
        last_prob = 0.0

        for _, row in sub.iterrows():
            pos = int(row["pos"])
            vals = [float(row[s]) for s in STATES]
            idx = max(range(3), key=lambda i: vals[i])
            state, prob = STATES[idx], vals[idx]
            if last_state is None:
                last_state, last_prob, last_pos = state, prob, pos
                continue

            gap = pos - last_pos
            if gap > max_gap:
                rows.append([chrom, start, last_pos, last_state if last_prob >= threshold else "NA"])
                rows.append([chrom, last_pos, pos, "NA"])
                start = pos
                last_state, last_prob, last_pos = state, prob, pos
                continue

            if last_prob >= threshold and prob >= threshold:
                if state == last_state:
                    last_pos = pos
                else:
                    mid = last_pos + round((pos - last_pos) / 2)
                    rows.append([chrom, start, mid, last_state])
                    start = mid
                    last_state, last_prob, last_pos = state, prob, pos
            elif last_prob < threshold and prob >= threshold:
                rows.append([chrom, start, pos, "NA"])
                start = pos
                last_state, last_prob, last_pos = state, prob, pos
            elif last_prob >= threshold and prob < threshold:
                rows.append([chrom, start, last_pos, last_state])
                start = last_pos
                last_state, last_prob, last_pos = state, prob, pos
            else:
                last_pos, last_state, last_prob = pos, state, prob

        rows.append([chrom, start, chrom_end, last_state if last_prob >= threshold else "NA"])

    merged = []
    for row in rows:
        if merged and merged[-1][0] == row[0] and merged[-1][3] == row[3] and merged[-1][2] >= row[1]:
            merged[-1][2] = max(merged[-1][2], row[2])
        else:
            merged.append(row)
    return merged


def main(config_path):
    c = configparser.ConfigParser(); c.read(config_path)
    root = Path(config_path).resolve().parent
    indir = root / c["OUTPUT"]["PREPARED_DIR"]
    outdir = root / c["OUTPUT"]["HMM_DIR"]
    outdir.mkdir(parents=True, exist_ok=True)
    lengths = load_chromosomes(root / c["INPUT"]["CHROM_FILE"])
    threshold = float(c["PARAMETER"]["PROBABILITY_THRESHOLD"])
    max_gap = int(c["PARAMETER"]["MAX_GAP"])

    for f in sorted(indir.glob("*.posterior")):
        df = pd.read_csv(f, sep="\t", dtype={"chrom": str})
        blocks = convert(df, lengths, threshold, max_gap)
        sample = f.name.rsplit(".posterior", 1)[0]
        out = outdir / f"{sample}_ahmm2bed"
        with open(out, "w") as g:
            g.write("chrom\tchromStart\tchromEnd\tancestry\n")
            for row in blocks:
                g.write("\t".join(map(str, row)) + "\n")
        print(f"Wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python 02.posterior_to_ahmm.py config.ini")
    main(sys.argv[1])
