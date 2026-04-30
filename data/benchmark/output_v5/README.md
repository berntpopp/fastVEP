# v5 — Current run (full SA stack + bug fixes)

This is the **current run** with all SA sources loaded and the SVI
wiring + indel-allele bugs fixed.

## What was loaded

| SA source | Loaded? | Notes |
|-----------|:-------:|-------|
| ClinVar (.osa)            | ✅ | 4,402,501 records |
| ClinVar protein (.oga)    | ✅ | 4,554 genes |
| gnomAD v4.1 exomes (.osa, per-chrom)         | ✅ | 25 chromosomes |
| gnomAD v4.1 gene constraints (.oga)         | ✅ | 18,173 genes |
| REVEL v1.3 (.osa, per-chrom)        | ✅ | 24 chromosomes |
| **PhyloP** (.osa, per-chrom)        | ✅ | distilled from gnomAD v4 INFO `phylop` (Zoonomia 241-mammal score) |
| **SpliceAI** (.osa, per-chrom)      | ✅ | distilled from gnomAD v4 INFO `spliceai_ds_max` |
| **ClinGen Gene-Disease Validity (.oga)**     | ✅ | 2,419 Definitive/Strong/Moderate genes — preferred over OMIM per ClinGen SVI / Abou Tayoun 2018 |

## Code fixes vs v1

1. **SpliceAI camelCase mismatch** in `sa_extract.rs`: writer wrote
   `spliceAI` but classifier matched `spliceai|spliceAi|splice_ai`.
   Added `spliceAI` to the match arms.
2. **PhyloP routing** in `sa_extract.rs`: PhyloP was attached to
   `aa.supplementary` (allele-level) but the classifier read it from
   `variant_supplementary` (variant-level). Read from both.
3. **Indel allele matching** in
   `analysis/acmg_benchmark/real_data/03_evaluate_concordance.py`:
   added `vep_allele(ref, alt)` to convert VCF (REF, ALT) →
   VEP CSQ Allele convention (`-` for deletion; insertion = right
   portion only) before grouping CSQ entries by allele. Without this
   all 48,539 indels in the truth fell into NoCall.

## Headline metrics

| Metric | Value | Δ vs v1 |
|--------|------:|--------:|
| Same-direction concordance | **65.1 %** | +10.4 pp |
| Exact match | **56.0 %** | +3.3 pp |
| Opposite direction | 0.06 % | (unchanged ~0) |
| Likely_benign recall | **42.4 %** | **+39 pp** |
| Benign recall | **58.0 %** | +25 pp |
| Pathogenic recall | **20.6 %** | +5 pp |
| Likely_pathogenic recall | **26.7 %** | +6 pp |

## Files

- `clinvar_2star.fastvep.vcf.gz` — bgzipped VCF with ACMG in CSQ INFO
  (~70 MB; vs ~25 GB for the prior pretty-printed JSON output)
- `concordance_matrix.csv` — 5-class truth × predicted matrix
- `concordance_summary.txt` — text rollup
- `concordance_by_chrom.csv` — per-chromosome breakdown
- `concordance_by_consequence.csv` — top consequences × class
- `criterion_firing_rates.csv` — per-criterion fire counts by truth class
- `rule_distribution.csv` — top criteria-set signatures
- `discrepancies.tsv` — opposite-direction calls (top 10k)
- `figures/` — 6 PNG + PDF figures (incl. v1 vs v5 comparison panels)

For the v1 baseline (no PhyloP / no SpliceAI / no ClinGen GDV), see
`../output_v1/`. For the version-by-version SA stack and code-fix
diff, see `../RUN_VERSIONS.md`.
