#!/usr/bin/env python3
"""Attach AHMM ancestry-block information to posterior sites."""

import configparser
import sys
from pathlib import Path


def main(config_path):
    c = configparser.ConfigParser(); c.read(config_path)
    root = Path(config_path).resolve().parent
    pdir = root / c["OUTPUT"]["PREPARED_DIR"]
    adir = root / c["OUTPUT"]["HMM_DIR"]
    outdir = root / c["OUTPUT"]["MERGED_DIR"]
    outdir.mkdir(parents=True, exist_ok=True)

    for pf in sorted(pdir.glob("*.posterior")):
        sample = pf.name.rsplit(".posterior", 1)[0]
        af = adir / f"{sample}_ahmm2bed"
        if not af.exists():
            raise FileNotFoundError(af)

        blocks = {}
        with open(af) as f:
            next(f, None)
            for line in f:
                p = line.rstrip().split("\t")
                if len(p) >= 4:
                    blocks.setdefault(p[0], []).append((int(p[1]), int(p[2]), p[3], f"{p[0]}-{p[1]}-{p[2]}", int(p[2]) - int(p[1])))

        out = outdir / f"{sample}.merged.posterior.bed.txt"
        with open(pf) as fin, open(out, "w") as fout:
            header = fin.readline().strip()
            fout.write(header + "\tancestry\tid\tlength\n")
            for line in fin:
                p = line.rstrip().split("\t")
                chrom, pos = p[0], int(p[1])
                anc = bid = length = "NA"
                for s, e, a, bi, l in blocks.get(chrom, []):
                    if s <= pos <= e:
                        anc, bid, length = a, bi, str(l)
                        break
                fout.write(line.rstrip() + f"\t{anc}\t{bid}\t{length}\n")
        print(f"Wrote {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python 03.merge_ancestry_info.py config.ini")
    main(sys.argv[1])
