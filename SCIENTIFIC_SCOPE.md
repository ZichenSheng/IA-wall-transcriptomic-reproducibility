# Scientific scope

## Evidence across cohorts

The literature audit covers 451 publication families and 1,345 HGNC genes. Figure 2 evaluates 239 gene–endpoint claims against genes matched for expression and variability in cohorts unused by the source papers. The matched percentile measures collective directional information. No individual claim in the 239-claim analysis survived Benjamini–Hochberg correction. The separate strict set has 58 publication-specific claims involving 49 genes. Its 30 overlapping rows map to 29 distinct collective units because two APP claims map to the same unit. Nominal calibration in the strict set follows its own definitions.

## Tissue compartments and interpretation

The lesion–rupture analysis uses 88 specimens across four rupture cohorts. Its six pairwise correlations share cohorts. The GSE54083 sensitivity removes pairs involving its rupture contrast while retaining the corrected unruptured aneurysm versus superficial temporal artery contrast in the lesion reference. The archived one-sided P value for this sensitivity is approximately 0.055. The platelet-, neutrophil- and erythroid-associated signal describes tissue sampled after haemorrhage. Pre-rupture attribution requires prospectively timed measurements. Wall-programme estimates depend on the control artery.

The nine published scores frequently reflect shared lesion, blood, clot, myeloid, interferon or brain-associated RNA. Performance in independent tissue and specificity for a named mechanism require separate evidence. Two scores lack a required high-weight feature for projection. The reproducibility matrix retains evaluations with limited power, unavailable inputs or failed calibration.

The human atlas and Xenium analyses use donors as inferential units and cover unruptured aneurysms. Attenuated contractile RNA and unstable transferred labels support compartment-aware interpretation. Composition, state and per-cell RNA output can all contribute to these patterns. Spatial centroid proximity describes the arrangement of cells at the measured resolution.

## What the simulations identify

The manuscript uses virtual tissues based on 21 unruptured aneurysm RNA-sequencing specimens to distinguish detection from biological attribution. Contribution scaling and uniform changes in per-cell RNA output share the same operator in the model. State substitution can yield a similar deconvolved smooth muscle cell (SMC) RNA share, and switching control references can shift that share even without a disease state substitution. Matching in the noiseless model equates mean deconvolution readings by construction. The reading estimates an RNA mixture coefficient, so absolute cell number, state and per-cell output remain unresolved in individual tissue.

## What the software runs

From prepared inputs, the release recomputes the collective analysis and its matched null and 31-cluster bootstrap, the observed lesion–rupture statistic and its GSE54083 sensitivity, and the correlation and sign agreement of paired GSE54083 gene effects. The GSE54083 comparison starts downstream of probe reconstruction and effect estimation. The joint phenotype-label permutation and rupture simulation tails use archived draws. Selected additional checks use frozen outputs, including summaries from the atlas and virtual-tissue analyses. The full simulation inputs and upstream raw-data workflows are outside the release. [Input provenance](DATA_SOURCES.md) · [Analysis-by-analysis coverage](REPRODUCIBILITY_MATRIX.tsv)
