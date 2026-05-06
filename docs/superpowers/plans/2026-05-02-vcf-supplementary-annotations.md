# VCF Supplementary Annotation Projection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Emit VCF-compatible non-JSON INFO projections for fastVEP's official supplementary annotation sources.

**Architecture:** Add a small VCF projection layer to `fastvep-io` that owns INFO headers, INFO replacement, escaping, and source-specific projection. Keep `pipeline.rs` focused on orchestration by passing loaded source keys and annotated variants to the projection layer.

**Tech Stack:** Rust 2021, `serde_json`, existing `fastvep-io`, `fastvep-cli`, and `fastvep-sa` crates.

---

### Task 1: Add Failing Projection Unit Tests

**Files:**
- Modify: `crates/fastvep-io/src/output.rs`

- [x] Add tests that build a small `VariationFeature` with supplementary keys for ClinVar, gnomAD, dbSNP, COSMIC, 1000G, TOPMed, MitoMap, PhyloP, GERP, DANN, REVEL, PrimateAI, dbNSFP, SpliceAI, OMIM, gnomAD gene constraints, and ClinVar protein.
- [x] Assert the projected INFO fields include `FV_*` fields and standard `SpliceAI`.
- [x] Assert projected values contain no raw `{`, `}`, or unescaped JSON quotes.
- [x] Run:

```bash
PATH="$HOME/.cargo/bin:$PATH" cargo test -p fastvep-io vcf_projection -- --nocapture
```

Expected: tests fail because the projection API does not exist yet.

### Task 2: Implement Projection Layer

**Files:**
- Modify: `crates/fastvep-io/src/output.rs`

- [x] Add projection definitions for all official fastSA source keys.
- [x] Add shared VCF subfield escaping.
- [x] Add `format_supplementary_vcf_info(&VariationFeature) -> Vec<(String, String)>`.
- [x] Keep `format_spliceai_info` behavior but route it through the projection layer.
- [x] Run the `fastvep-io` tests and confirm Task 1 passes.

### Task 3: Add Failing Header and INFO Replacement Tests

**Files:**
- Modify: `crates/fastvep-cli/tests/sa_build_oga.rs`

- [x] Add an end-to-end test where the input VCF already has `CSQ` and `SpliceAI` headers/values.
- [x] Assert output has one `CSQ` header, one `SpliceAI` header, and current fastVEP values replace stale input values.
- [x] Run:

```bash
PATH="$HOME/.cargo/bin:$PATH" cargo test -p fastvep-cli annotate_vcf_replaces_existing_fastvep_info -- --nocapture
```

Expected: test fails because `CSQ` is duplicated today.

### Task 4: Move Header and INFO Merge Policy Out of Pipeline

**Files:**
- Modify: `crates/fastvep-io/src/output.rs`
- Modify: `crates/fastvep-cli/src/pipeline.rs`

- [x] Add helpers for `##INFO` ID extraction and header replacement.
- [x] Add a helper that formats final INFO by replacing fastVEP-owned IDs and appending new projections.
- [x] Update `pipeline.rs` to request headers from `fastvep-io` using loaded SA and gene source keys.
- [x] Update `write_vcf_line` to call the shared INFO formatter.
- [x] Run the CLI replacement test and existing SpliceAI test.

### Task 5: Fix Source JSON Construction Risks

**Files:**
- Modify: `crates/fastvep-sa/src/sources/spliceai.rs`
- Optionally modify unsafe source builders discovered during implementation.

- [x] Add a failing SpliceAI parser test with a gene symbol containing VCF-special characters.
- [x] Use `serde_json::json!` or equivalent safe serialization for SpliceAI annotation records.
- [x] Run:

```bash
PATH="$HOME/.cargo/bin:$PATH" cargo test -p fastvep-sa sources::spliceai -- --nocapture
```

Expected: parser tests pass and generated JSON parses with `serde_json`.

### Task 6: Update Docs

**Files:**
- Modify: `README.md`

- [x] Update VCF output documentation to state VCF includes `CSQ`, standard `SpliceAI` when loaded, and fastVEP-specific `FV_*` supplementary INFO fields.
- [x] State JSON remains the richest structured output.

### Task 7: Full Verification and Commit

**Files:**
- All modified files.

- [x] Run:

```bash
PATH="$HOME/.cargo/bin:$PATH" cargo test --workspace
git diff --check
```

- [x] Validate a generated sample VCF with `pysam.VariantFile`.
- [x] Validate a larger real SpliceAI SNV+indel subset from indexed masked hg38 source VCFs; confirm `SpliceAI` is generated for all 400 records from fastSA, not passed through from input.
- [x] Add an indel regression proving allele-matched fastSA uses uploaded VCF alleles while positional fastSA keeps normalized coordinates.
- [x] Commit and push:

```bash
git add docs/superpowers/specs/2026-05-02-vcf-supplementary-annotations-design.md \
  docs/superpowers/plans/2026-05-02-vcf-supplementary-annotations.md \
  README.md crates/fastvep-io/src/output.rs crates/fastvep-cli/src/pipeline.rs \
  crates/fastvep-cli/tests/sa_build_oga.rs crates/fastvep-cli/fixtures/spliceai/gnomad-mini.vcf \
  crates/fastvep-sa/src/sources/spliceai.rs
git commit -m "Emit fastSA annotations as VCF INFO projections"
git push origin stream-spliceai-sa-build
```
