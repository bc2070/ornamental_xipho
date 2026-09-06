# Unified LKhard / LKsite Pipeline

This directory contains a unified pipeline for the LKhard and LKsite methods.
The same scripts are used for both methods; method-specific behavior is selected
through `config.ini`.


## Select the method

Edit:

```ini
[GENERAL]
METHOD = LKhard
```

or:

```ini
[GENERAL]
METHOD = LKsite
```

Then submit the same SLURM script:

```bash
mkdir -p logs
sbatch main_pipeline.sh
```

## Pipeline structure

```text
00.prepare_posterior.py
        |
        v
01.calculate_likelihood.py
        |
        v
02.posterior_to_ahmm.py
        |
        v
03.merge_ancestry_info.py
        |
        v
04.randomization.py
        |
        v
01.calculate_likelihood.py
        |
        v
05.baseline.py
        |
        v
06.plot_manhattan.R
```

`00.prepare_posterior.py` is shared by both methods. It detects two- or
three-ancestry posterior files, ranks ancestries by genome-wide contribution,
and for three ancestries retains the two largest contributions and maps them
to the standard states `2,0`, `1,1`, and `0,2`.

`01.calculate_likelihood.py` contains the method-specific likelihood formulas.
For LKhard it uses HWE genotype probabilities; for LKsite it uses Eq. S1 and
the minor-parental contribution `h`.

`04.randomization.py` keeps the method-specific null models: tract placement
for LKhard and chromosome-block permutation for LKsite.

All comments and documentation inside scripts are written in English.

## Method-specific behavior

### LKhard

- Uses HWE genotype probabilities derived from the standardized posterior.
- Converts posterior probabilities into ancestry blocks using the common AHMM conversion step.
- Randomizes ancestry tracts while preserving tract lengths.
- Calculates the null distribution from randomized tract configurations.

### LKsite

- Estimates `h` as the genome-wide contribution of the minor parental ancestry.
- Uses the LKsite Eq. S1 site-level likelihood.
- Randomizes chromosome ancestry blocks by permutation while retaining the original block structure.
- Calculates the null distribution from the randomized site-level scores.

The two methods therefore share the same input preparation, block conversion,
merging, threshold calculation, and plotting framework, while retaining their
original method-specific likelihood and randomization models.
