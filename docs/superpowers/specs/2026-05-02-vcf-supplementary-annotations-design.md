# VCF Supplementary Annotation Projection Design

## Goal

Make `fastvep annotate --output-format vcf` emit VCF-compatible annotations for fastVEP's official supplementary annotation sources without writing JSON payloads into the INFO column.

## Context

fastSA annotation providers expose annotations to the annotation pipeline as `json_key` plus pre-serialized JSON strings. That storage/query boundary is intentional for JSON output and classification, and `.osa2` may store typed arrays internally before reconstructing JSON at the provider boundary.

VCF output is different. VCF INFO fields require declared headers, declared cardinality and type, delimiter-safe values, and source-specific projection. The existing writer emits VEP-style `CSQ` and now emits standard SpliceAI `SpliceAI`, but other fastSA sources are available in JSON output only.

## Design

Add a reusable VCF projection layer in `fastvep-io` that owns:

- INFO header definitions.
- Replacement of existing fastVEP-owned INFO headers.
- Replacement of existing fastVEP-owned INFO values.
- Shared delimiter escaping for VCF INFO strings.
- Source-specific conversion from fastSA JSON values to pipe-delimited INFO strings.

The CLI pipeline should only assemble available sources and call this layer. It should not know per-source header or INFO formatting details.

## Output Schema

Keep established external conventions where they exist:

- `CSQ`: VEP-style consequence annotation.
- `SpliceAI`: standard SpliceAI format:
  `ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`.

For fastVEP-specific supplementary projections, use declared `FV_*` INFO fields. Each field is `Number=.` and `Type=String`, with comma-separated entries and pipe-delimited subfields. Examples:

- `FV_CLINVAR=G|Pathogenic|criteria_provided|Breast_cancer|SNV|SO%3A0001483`
- `FV_GNOMAD=G|1.200000e-04|12|100000|0|...`
- `FV_REVEL=G|0.8123`
- `FV_OMIM=BRCA1|113705|Breast%20cancer&ovarian%20cancer`

Every emitted field must have a matching `##INFO=<ID=...>` line. Values must not contain raw JSON braces, quotes, semicolons, tabs, newlines, unescaped spaces, or unescaped commas inside a subfield.

For VCF-backed allele-specific fastSA sources, lookups and emitted `ALLELE` subfields use the uploaded VCF POS/REF/ALT, not fastVEP's minimized internal allele representation. This is required for indels, where the internal allele may be `-` or a stripped insertion sequence while the source `.osa` is keyed by the original VCF record. Positional, non-allele-matched sources continue to query fastVEP's normalized variant coordinate.

## Source Coverage

Variant/allele-level `.osa` sources:

- ClinVar -> `FV_CLINVAR`
- gnomAD -> `FV_GNOMAD`
- dbSNP -> `FV_DBSNP`
- COSMIC -> `FV_COSMIC`
- 1000 Genomes -> `FV_1KG`
- TOPMed -> `FV_TOPMED`
- MitoMap -> `FV_MITOMAP`
- PhyloP -> `FV_PHYLOP`
- GERP -> `FV_GERP`
- DANN -> `FV_DANN`
- REVEL -> `FV_REVEL`
- SpliceAI -> `SpliceAI`
- PrimateAI -> `FV_PRIMATEAI`
- dbNSFP -> `FV_DBNSFP`

Gene-level `.oga` sources:

- OMIM -> `FV_OMIM`
- gnomAD gene constraints -> `FV_GNOMAD_GENE`
- ClinVar protein -> `FV_CLINVAR_PROTEIN`

## Replacement Policy

fastVEP owns these emitted IDs during annotation. If an input VCF already contains one of these headers or INFO values and fastVEP is going to emit a value for the same ID, fastVEP replaces the old one. This follows Ensembl VEP's default behavior for `CSQ`.

If a source is not loaded or a variant has no value for a source, no source-specific INFO value is emitted. Headers should only be emitted for `CSQ` plus loaded/available projection sources.

## Testing

Tests must cover:

- `CSQ` header and INFO replacement when input already has `CSQ`.
- Existing `SpliceAI` replacement and no duplicate header.
- VCF projection for all supported fastSA keys using unit-level `VariationFeature` fixtures.
- No raw JSON delimiters in projected `FV_*` INFO values.
- Special-character escaping in subfields.
- SpliceAI source JSON construction with special characters in gene symbols.
- End-to-end CLI output remains parseable by a VCF parser for representative fixtures.
- Indel fixtures prove allele-matched SpliceAI lookup and VCF `ALLELE` output use uploaded VCF alleles while positional scores still use normalized indel coordinates.

## Non-Goals

- Do not dump arbitrary JSON into VCF.
- Do not claim SnpEff `ANN` compatibility in this change.
- Do not invent support for sources not already supported by `sa-build`.
- Do not change JSON output structure except where invalid source JSON construction is corrected.
