# Reproducible molecular information and the limits of biological identifiability in intracranial aneurysm wall transcriptomics

**Zichen Sheng, Mingze Wei, Kunjian Lei and Xingen Zhu**

## Study overview

Intracranial aneurysm (IA) studies report many molecular markers and named pathways, often drawing on the same cohorts. This study asks which findings reproduce in tissue from cohorts their source studies had not used, and which biological processes the RNA measurements can identify. Published gene directions reproduced collectively and mapped to tissue programmes. Aneurysm walls resected after rupture showed a haemostatic and leukocyte-associated signal beyond the shared lesion axis. Published scores, human single-cell and spatial data, and simulations with known generating processes showed why these signals can be informative about tissue compartments while leaving their cellular origins and mechanisms uncertain. The findings motivate anatomically matched controls, measurements that resolve tissue compartments, and prospective sampling for studies of rupture risk.

## Main contributions

- **Collective evidence from published claims.** An audit of 451 publication families identified 239 evaluable gene–endpoint claims. In cohorts unused by their source studies, the reported directions achieved a mean expression-matched percentile of 0.738 against a null mean of 0.499. No individual claim survived multiple-testing correction.
- **Lesion and rupture-associated RNA.** A shared lesion axis captured reduced contractile RNA, matrix gain and myeloid enrichment. Walls sampled after haemorrhage showed additional platelet-, neutrophil- and erythroid-associated RNA in proportions distinct from circulating blood. The choice of control artery shaped estimates of wall programmes. Establishing which signals precede rupture requires prospective evidence.
- **What transcriptomic scores measure.** Nine published multigene scores often tracked shared lesion, blood, clot, immune or brain-associated RNA. Their ability to distinguish tissue groups therefore had limited value as evidence for their named mechanisms.
- **Biological and feature identity.** Human single-cell and spatial data showed reduced contractile RNA across the mural-like population and instability in transferred labels. Simulations showed that smooth muscle cell state changes, reduced RNA contribution and reduced RNA output can produce overlapping deconvolution readings. At the feature level, reconstructing probe identities in GSE54083 revealed that the default and corrected gene-level effects were nearly uncorrelated across 18,918 genes.

## Reproducibility

The executable workflows use bundled, documented gene-level inputs to recompute the collective matched-percentile analysis and cluster bootstrap, the observed lesion–rupture excess correlations, and the GSE54083 annotation comparison. The runner also checks selected results against archived outputs.

Use Python 3.13 and run these commands from the repository root.

```sh
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m scripts.verify_source_data
python -m unittest discover -s tests -v
```

## Citation

Cite the companion manuscript, *Reproducible molecular information and the limits of biological identifiability in intracranial aneurysm wall transcriptomics*, by Zichen Sheng, Mingze Wei, Kunjian Lei and Xingen Zhu. Cite the software version used with its version-specific DOI. [Software citation](CITATION.cff) · [GitHub v1.0.1](https://github.com/ZichenSheng/IA-wall-transcriptomic-reproducibility/releases/tag/v1.0.1) · [Zenodo v1.0.1](https://doi.org/10.5281/zenodo.24070166)

## Data sources and documentation

[Data sources and prepared inputs](DATA_SOURCES.md) · [Scientific scope](SCIENTIFIC_SCOPE.md) · [Reproducibility matrix](REPRODUCIBILITY_MATRIX.tsv) · [Environment](ENVIRONMENT.md) · [Provenance](PROVENANCE.md)
