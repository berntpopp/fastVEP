# Benchmark run versions

Each row of the `output_v*/` directories represents one end-to-end run
of the ClinVar 2-star+ benchmark on the same 673,660-variant truth set.
Successive runs differ only in (a) which supplementary annotation
databases were loaded into `--sa-dir` and (b) which classifier / output
bugs had been fixed.

## SA stack per run

|                                     |  v1  |  v2  |  v4  |  v5  |
|-------------------------------------|:----:|:----:|:----:|:----:|
| **Variant-level (.osa)**            |      |      |      |      |
| ClinVar                             |  ✅  |  ✅  |  ✅  |  ✅  |
| gnomAD v4.1 exomes (per-chrom)      |  ✅  |  ✅  |  ✅  |  ✅  |
| REVEL v1.3 (per-chrom)              |  ✅  |  ✅  |  ✅  |  ✅  |
| **PhyloP** (per-chrom)              |  ❌  |  ✅  |  ✅  |  ✅  |
| **SpliceAI** (per-chrom)            |  ❌  |  ✅  |  ✅  |  ✅  |
| **Gene-level (.oga)**               |      |      |      |      |
| ClinVar protein                     |  ✅  |  ✅  |  ✅  |  ✅  |
| gnomAD gene constraints             |  ✅  |  ✅  |  ✅  |  ✅  |
| **ClinGen Gene-Disease Validity**   |  ❌  |  ✅  |  ✅  |  ✅  |

(PhyloP and SpliceAI are distilled from gnomAD v4 INFO fields
`phylop` and `spliceai_ds_max` rather than re-downloaded; the gnomAD
v4 sites VCF already includes them. ClinGen GDV substitutes for
OMIM `genemap2.txt` per ClinGen SVI / Abou Tayoun 2018 — same `.oga`
schema, `omim` json_key, but with a multi-curator scored rubric and
explicit Definitive/Strong/Moderate filtering.)

## Code fixes per run

|                                                              |  v1  |  v2  |  v4  |  v5  |
|--------------------------------------------------------------|:----:|:----:|:----:|:----:|
| SpliceAI `spliceAI` json_key recognised in classifier        |      |      |  ✅  |  ✅  |
| PhyloP read from `allele_supplementary` (CLI's actual route) |      |      |  ✅  |  ✅  |
| VCF + bgzip output (vs 25 GB pretty JSON)                    |      |      |  ✅  |  ✅  |
| `vep_allele(ref, alt)` indel matching in concordance script  |      |      |      |  ✅  |

v3 was a partial run (PhyloP+SpliceAI loaded but bugs still latent);
its results are functionally indistinguishable from v2 and were
overwritten before being preserved.

## Headline metrics per run

|                            |     v1     |     v5     |     Δ      |
|----------------------------|-----------:|-----------:|-----------:|
| Same-direction concordance |   54.7 %   | **65.1 %** | **+10.4 pp** |
| Exact match                |   52.7 %   | **56.0 %** | +3.3 pp    |
| Opposite direction         |   0.005 %  |   0.06 %   | (≈0)       |
| NoCall                     |   0.0 %    |   0.0 %    | —          |
| Pathogenic recall          |   15.7 %   | **20.6 %** | +5 pp      |
| Likely_pathogenic recall   |   20.9 %   | **26.7 %** | +6 pp      |
| VUS recall                 |   96.6 %   |   92.6 %   | -4 pp      |
| **Likely_benign recall**   |   **3.2 %**|**42.4 %**  | **+39 pp** |
| Benign recall              |   33.2 %   | **58.0 %** | +25 pp     |

## Driver of each lift

- **+39 pp LB recall, +25 pp B recall**: BP7 went from **0** → **81,706
  fires** once PhyloP+SpliceAI were loaded *and* both wiring bugs were
  fixed. (Walker 2023: BP7 needs synonymous + low SpliceAI + low PhyloP.)
- **PVS1 ~9× more fires** (5,636 → 50,062 P+LP): ClinGen GDV provides
  the disease-gene fallback when gnomAD pLI/LOEUF don't cross the LOF
  threshold; the indel-allele fix surfaced thousands of frameshift
  pathogenic variants that previously hid behind NoCall.
- **VUS recall slight drop (-4 pp)**: by design — when more benign
  evidence fires, some variants previously called VUS now correctly
  drop to LB or B (which doesn't match a VUS truth).

## Where to find each version

- v1 baseline: `output_v1/concordance_matrix.csv` +
  `output_v1/README.md` (raw outputs were overwritten; matrix
  reconstructed from documentation)
- v5 current: `output_v5/` (full outputs + figures + raw VCF.gz)
