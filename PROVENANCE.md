# Provenance

This package supports the November 9, 2026 final IA manuscript, its Figure 1 V4 and the November 7 supervisor-corrected Figures 2–6. The final Supplementary Information and Tables S1–S5 supply the accompanying definitions and result tables. `expected/source_authority_sha256.tsv` records the final files and frozen analysis authorities. `expected/input_sha256.json` identifies the distributed inputs.

| Script | Frozen authority and adaptation |
|---|---|
| `_q1_kernel.py` | Oct 27 `p0_q1_core.py` and `p0_dependency.py`. Extracted DL pooling, percentile, replacement null and cluster bootstrap. |
| `reproduce_gene_claims.py` | Expression matching follows `Q1.expr_pool` / `signed_pool` in `p0_q1_core.py`. Frozen case definition and cohort eligibility are inputs. SciPy BH was checked against all archived claim q values. |
| `_d1_kernel.py` | Oct 11 decisive execution `d1lib.py` and `d1_engine.py`. Extracted one-factor fit, inverse-normal ranking, correlations and excess statistic. |
| `reproduce_rupture_decomposition.py` | Same D1 contrast order and `S4_DROP_RUPTURE_GSE54083` rule from `d1_sensitivity.py`. IO uses the prepared eight-vector projection. Archived inference checks follow `d1_primary.py` and `d1_sensitivity.py`. |
| `verify_source_data.py` | Release entry point, Figure 1E paired-effect calculation and frozen aggregate checks. |

The audit traced upstream cohort parsing, frozen annotation, effect estimation and variance, matching and DL pooling, permutations and bootstrap, D1, score definitions/transport, atlas summaries, and E01/M01 through manuscript Source Data. Score refitting, atlas/spatial rebuilding, wet-lab validation, the full joint label permutation, full M1 generation, the complete Fig. 6 simulation and full figure rendering are outside this archive. Fig. 6 uses the approved V9 numerical lineage with November 7 labels. Figure 1 V4 preserves the existing numerical lineage. Its paired-effect calculation uses the exact source table for the approved panel.

SHA-256 records identity and detects changes to files. Scholarly attribution, computational validation, biological interpretation and priority are distinct questions.
