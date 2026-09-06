#!/usr/bin/env python3
"""Generate randomized data using the method-specific null model."""

import configparser
import random
import sys
from collections import defaultdict
from pathlib import Path


def natural_key(x):
    return (0, int(x)) if str(x).isdigit() else (1, str(x))


def read_chrom_lengths(path):
    out = []
    with open(path) as f:
        for line in f:
            p = line.strip().split()
            if not p or p[0].lower() in {"chrom", "chr", "chromosome"}:
                continue
            out.append((p[0], int(float(p[2] if len(p) >= 3 else p[1]))))
    return out


def random_tracts(chroms, fragments):
    result = []
    for length, state in fragments:
        eligible = [(c, end) for c, end in chroms if end >= length]
        if not eligible:
            continue
        chrom, max_end = random.choice(eligible)
        start = random.randint(0, max_end - length)
        result.append((chrom, start, start + length, state))
    return result


def choose_ancestry(state_string):
    """Resolve overlapping ancestry states using the original priority rule."""
    states = set(state_string.split(","))
    priority = ["0,2", "1,1", "2,0"]
    for state in priority:
        if state in states:
            return state
    return "NA"


def merge_intervals(rows):
    rows = sorted(rows, key=lambda x: (natural_key(x[0]), x[1], x[2]))
    merged = []
    for chrom, start, end, state in rows:
        if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
            merged[-1][2] = max(merged[-1][2], end)
            merged[-1][3] = choose_ancestry(merged[-1][3] + "," + state)
        else:
            merged.append([chrom, start, end, state])
    return merged


def randomize_lkhard(c, root):
    n = int(c["PARAMETER"]["ITERATIONS"])
    chroms = read_chrom_lengths(root / c["INPUT"]["CHROM_FILE"])
    indir = root / c["OUTPUT"]["HMM_DIR"]
    outdir = root / c["OUTPUT"]["RANDOM_DIR"]
    outdir.mkdir(parents=True, exist_ok=True)

    for f in sorted(indir.glob("*_ahmm2bed")):
        sample = f.name.rsplit("_ahmm2bed", 1)[0]
        rows = []
        with open(f) as h:
            next(h, None)
            for line in h:
                p = line.rstrip().split("\t")
                if len(p) >= 4:
                    rows.append((p[0], int(p[1]), int(p[2]), p[3]))
        het = [(e - s, st) for _, s, e, st in rows if st == "1,1"]
        minor = [(e - s, st) for _, s, e, st in rows if st == "0,2"]
        background = [(ch, s, e, st) for ch, s, e, st in rows if st not in {"1,1", "0,2"}]

        for i in range(1, n + 1):
            randomized = background + random_tracts(chroms, het) + random_tracts(chroms, minor)
            randomized = merge_intervals(randomized)
            out = outdir / f"{sample}.rand.{i}.fixed.bed"
            with open(out, "w") as g:
                for row in randomized:
                    g.write("\t".join(map(str, row)) + "\n")


def randomize_lksite(c, root):
    n = int(c["PARAMETER"]["ITERATIONS"])
    indir = root / c["OUTPUT"]["MERGED_DIR"]
    outdir = root / c["OUTPUT"]["RANDOM_DIR"]
    outdir.mkdir(parents=True, exist_ok=True)

    for mf in sorted(indir.glob("*.merged.posterior.bed.txt")):
        sample = mf.name.rsplit(".merged.posterior.bed.txt", 1)[0]
        bychrom = defaultdict(list)
        with open(mf) as f:
            next(f, None)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 8:
                    bychrom[p[0]].append(p)

        blocks = {}
        for chrom, rows in bychrom.items():
            groups = defaultdict(list)
            for row in rows:
                groups[row[6]].append(row)
            blocks[chrom] = list(groups.values())

        for k in range(1, n + 1):
            allrows = []
            for chrom in sorted(blocks, key=natural_key):
                shuffled = blocks[chrom][:]
                random.shuffle(shuffled)
                source_rows = [r for block in shuffled for r in block]
                coords = [(r[0], r[1]) for r in bychrom[chrom]]
                info = [r[2:5] for r in source_rows]
                allrows.extend([list(coord) + x for coord, x in zip(coords, info)])

            out = outdir / f"{sample}.rand.all.{k}.txt"
            with open(out, "w") as g:
                g.write("chrom\tpos\t2,0\t1,1\t0,2\n")
                for row in allrows:
                    g.write("\t".join(row) + "\n")


def main(config_path):
    c = configparser.ConfigParser(); c.read(config_path)
    root = Path(config_path).resolve().parent
    method = c["GENERAL"]["METHOD"].strip()
    if method == "LKhard":
        randomize_lkhard(c, root)
    elif method == "LKsite":
        randomize_lksite(c, root)
    else:
        raise ValueError("METHOD must be LKhard or LKsite")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python 04.randomization.py config.ini")
    main(sys.argv[1])
