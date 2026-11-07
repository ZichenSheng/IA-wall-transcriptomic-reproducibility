# Provenance

This candidate corresponds to the clean manuscript and revised Supplementary Methods dated 2026-11-07, Supplementary Tables S1–S5, the six supervisor-revised PDFs and the corresponding frozen Source Data. See `expected/source_authority_sha256.tsv` for byte-level links to the authorities and `expected/input_sha256.json` for the distributed projections. The version is `1.0.0-rc.1`. The Git tag and external archive receipt identify the exact commit. The version DOI and article or preprint DOI are pending and belong in the final citation and archive receipt once verified.

| Script | Frozen authority and adaptation |
|---|---|
| `_q1_kernel.py` | Oct 27 `p0_q1_core.py` and `p0_dependency.py`. Extracted DL pooling, percentile, replacement null and cluster bootstrap. |
| `reproduce_gene_claims.py` | Expression matching follows `Q1.expr_pool` / `signed_pool` in `p0_q1_core.py`. Frozen case definition and cohort eligibility are inputs. SciPy BH was checked against all archived claim q values. |
| `_d1_kernel.py` | Oct 11 decisive execution `d1lib.py` and `d1_engine.py`. Extracted one-factor fit, inverse-normal ranking, correlations and excess statistic. |
| `reproduce_rupture_decomposition.py` | Same D1 contrast order and `S4_DROP_RUPTURE_GSE54083` rule from `d1_sensitivity.py`. IO uses the prepared eight-vector projection. Archived inference checks follow `d1_primary.py` and `d1_sensitivity.py`. |
| `reproduce_smc_simulation.py` | Original E01 `e01_refs.py` and Fig. 6 M01 `m01_matched_estimate_residuals.py`. Selected full-dose worlds. The linear per-gene operator and library sum are represented on the sufficient NNLS inputs. Verified against the archived M01 values used by the final Fig. 6 lineage. |
| `verify_source_data.py` | Release entry point and Level C checks written for this candidate. Checks frozen aggregate values and links them to the implemented computations when applicable. |

The audit traced upstream cohort parsing, frozen annotation, effect estimation and variance, matching and DL pooling, permutations and bootstrap, D1, score definitions/transport, atlas summaries, and E01/M01 through manuscript Source Data. Score refitting, atlas/spatial rebuilding, wet-lab validation, the full joint label permutation, full M1 generation, the complete Fig. 6 simulation and full figure rendering are outside this archive. Fig. 6 uses the approved V9 numerical lineage with November 7 labels. The experimental figure redesigns are outside the authoritative lineage.

SHA-256 records identity and detects changes to files. Scholarly attribution, computational validation, biological interpretation and priority are distinct questions.
