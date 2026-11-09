# Reproducible molecular information and the limits of biological identifiability in intracranial aneurysm wall transcriptomics

Zichen Sheng, Mingze Wei, Kunjian Lei and Xingen Zhu

Companion software and prepared inputs for the final IA manuscript dated November 9, 2026. This package evaluates the reproducibility of published gene directions and the extent to which tissue RNA can identify underlying biology in human intracranial aneurysm (IA).

## Run

Use Python 3.13 from the repository root.

```sh
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m scripts.verify_source_data
python -m unittest discover -s tests -v
```

The runner returns JSON on success and exits nonzero on input or numerical mismatch. The standalone run needs only the package contents and installed Python dependencies. Use `--workflow claims` or `--workflow rupture` to run one analysis. See [environment](ENVIRONMENT.md) for the tested installation.

## What the package reproduces

- **Figure 1E, Level B.** Recomputes Pearson correlation and sign agreement from 18,918 paired, prepared GSE54083 Hedges’ g values under the default and reconstructed annotations. The expected results are r = −0.0082494 and 50.29% sign agreement.

- **Figure 2A/B, Level B.** Starting with gene-by-cohort Hedges’ g, sampling variances and expression/variability ranks, recomputes cohort pooling, 200 matched comparators per claim, all 239 matched percentiles, matched-null draws and the 31-cluster bootstrap. The 58 publication-specific strict claims are a separate Level C check.
- **Figure 3A, Level B.** Starting with eight prepared Welch vectors, refits the lesion factor, calculates six pairwise excess correlations and averages them, then removes rupture pairs involving GSE54083 for one sensitivity. GSE54083 uses the reconstructed IA annotation. The Bayesian intervals and Monte Carlo P values are checked from archived draws at Level C.
- **Other Level C checks.** Verifies selected Figure 1, 3B, 4, 5 and 6 aggregates and the archived Figure 2 joint label-permutation tail. See the [reproducibility matrix](REPRODUCIBILITY_MATRIX.tsv) for the result-by-result boundary.

Level A starts with public repository files. Level B starts with documented prepared inputs. Level C checks existing outputs. This package implements B and C. The [data sources](DATA_SOURCES.md) explain the input preparation and the [scope note](SCIENTIFIC_SCOPE.md) describes the interpretation limits.

## Citation and access

Use the ordered author list and version in [CITATION.cff](CITATION.cff). The [Zenodo version record](https://doi.org/10.5281/zenodo.24069701) identifies this archive. The original code and documentation use the MIT license. Data use remains subject to the applicable source terms. See [provenance](PROVENANCE.md).
