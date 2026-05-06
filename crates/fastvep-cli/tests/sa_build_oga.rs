//! End-to-end tests for supplementary annotation CLI paths.
//!
//! Each `sa-build` test writes a small fixture for the source, calls
//! `run_sa_build` (the same entrypoint the CLI uses), and reads the resulting
//! database back to confirm the round-trip.

use fastvep_cli::pipeline::{run_annotate, run_sa_build, AnnotateConfig};
use fastvep_sa::gene::GeneIndex;
use std::fs::{self, File};
use std::io::Write;

const SPLICEAI_SOURCE_VCF: &str = include_str!("../fixtures/spliceai/spliceai-mini.vcf");
const SPLICEAI_INDEL_SOURCE_VCF: &str =
    include_str!("../fixtures/spliceai/spliceai-indel-mini.vcf");
const GNOMAD_SOURCE_VCF: &str = include_str!("../fixtures/spliceai/gnomad-mini.vcf");
const INPUT_NO_SPLICEAI_INFO_VCF: &str =
    include_str!("../fixtures/spliceai/input-no-spliceai-info.vcf");
const MINI_GFF3: &str = include_str!("../fixtures/spliceai/mini.gff3");

fn read_oga(path: &std::path::Path) -> GeneIndex {
    let mut f = File::open(path).expect("open .oga");
    GeneIndex::read_from(&mut f).expect("parse .oga")
}

#[test]
fn sa_build_omim_writes_oga_with_records() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("genemap2.txt");
    let output = tmp.path().join("omim");

    // Minimal genemap2.txt fixture — column layout matches the real OMIM
    // file format (gene symbol at index 5, MIM at index 8, phenotypes at 12).
    let fixture = "# Generated\n\
                   # Copyright OMIM\n\
                   1\tp36.33\t1:10001-20000\tGene1\t\tBRCA1\tprotein\t\t113705\t\t\t\tBreast cancer, 114480 (3), Autosomal dominant; Ovarian cancer, 167000 (3)\n\
                   1\tp36.33\t1:30001-40000\tGene2\t\tTP53\tprotein\t\t191170\t\t\t\tLi-Fraumeni syndrome 1, 151623 (3), Autosomal dominant\n";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "omim",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    let oga_path = output.with_extension("oga");
    assert!(oga_path.exists(), ".oga file should be written");

    let index = read_oga(&oga_path);
    assert_eq!(index.header.json_key, "omim");
    assert_eq!(index.header.name, "ClinGen GDV / OMIM");
    assert!(index.gene_count() >= 2);

    let brca1 = index.get("BRCA1").expect("BRCA1 should be present");
    let json = brca1.first().unwrap();
    assert!(
        json.contains("113705"),
        "BRCA1 should have MIM 113705: {}",
        json
    );
    assert!(json.contains("Breast cancer"));
}

#[test]
fn sa_build_gnomad_genes_writes_oga_with_records() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("constraint.tsv");
    let output = tmp.path().join("gnomad_genes");

    let fixture = "\
gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tmis_z\tsyn_z
BRCA1\tENST00000357654\t0\t50.2\t0.00\t0.03\t1.0000\t3.45\t0.12
TP53\tENST00000269305\t0\t25.1\t0.00\t0.05\t0.9999\t5.67\t-0.34
";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "gnomad_genes",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    let oga_path = output.with_extension("oga");
    let index = read_oga(&oga_path);
    assert_eq!(index.header.json_key, "gnomad_genes");
    assert_eq!(index.gene_count(), 2);

    let brca1 = index.get("BRCA1").unwrap();
    let json = brca1.first().unwrap();
    assert!(json.contains("\"pLI\":1.0000"));
    assert!(json.contains("\"loeuf\":0.0300"));
    assert!(json.contains("\"misZ\":3.45"));
}

#[test]
fn sa_build_gnomad_gene_alias_routes_to_same_key() {
    // Both `gnomad_genes` (plural) and `gnomad_gene` (singular) are accepted
    // for the CLI; both must produce a database with json_key=gnomad_genes
    // so the classifier picks them up consistently.
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("constraint.tsv");
    let output = tmp.path().join("gnomad_gene_alias");
    let fixture = "\
gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tmis_z\tsyn_z
BRCA1\tENST00000357654\t0\t50.2\t0.00\t0.03\t1.0000\t3.45\t0.12
";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "gnomad_gene",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    let index = read_oga(&output.with_extension("oga"));
    assert_eq!(index.header.json_key, "gnomad_genes");
}

#[test]
fn sa_build_clinvar_protein_writes_oga_with_records() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("clinvar.vcf");
    let output = tmp.path().join("clinvar_protein");

    // Minimal ClinVar VCF: two pathogenic missense entries with protein
    // change in CLNHGVS, one rejected (Benign).
    let fixture = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t7676154\t12345\tG\tA\t.\t.\tCLNSIG=Pathogenic;MC=SO:0001583|missense_variant;GENEINFO=TP53:7157;CLNHGVS=NP_000537.3:p.Arg175His
17\t7676156\t12346\tT\tC\t.\t.\tCLNSIG=Likely_pathogenic;MC=SO:0001583|missense_variant;GENEINFO=TP53:7157;CLNHGVS=NP_000537.3:p.Arg248Trp
17\t7676160\t12347\tG\tA\t.\t.\tCLNSIG=Benign;MC=SO:0001583|missense_variant;GENEINFO=TP53:7157;CLNHGVS=NP_000537.3:p.Pro72Leu
";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "clinvar_protein",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    let index = read_oga(&output.with_extension("oga"));
    assert_eq!(index.header.json_key, "clinvar_protein");
    assert!(index.gene_count() >= 1);

    let tp53 = index.get("TP53").expect("TP53 should be present");
    let json = tp53.first().unwrap();
    // The two pathogenic entries should make it through; the Benign one shouldn't.
    assert!(
        json.contains("\"pos\":175"),
        "should include p.Arg175His: {}",
        json
    );
    assert!(
        json.contains("\"pos\":248"),
        "should include p.Arg248Trp: {}",
        json
    );
    assert!(
        !json.contains("\"pos\":72"),
        "Benign p.Pro72Leu must NOT be in index: {}",
        json
    );
}

#[test]
fn sa_build_unknown_source_errors_with_full_supported_list() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("nope.txt");
    File::create(&input).unwrap().write_all(b"").unwrap();

    let err = run_sa_build(
        "definitely_not_a_source",
        input.to_str().unwrap(),
        "out",
        "GRCh38",
    )
    .expect_err("must error on unknown source");
    let msg = format!("{}", err);
    // Error should list the new gene-level sources alongside variant-level ones.
    assert!(msg.contains("omim"), "error should mention omim: {}", msg);
    assert!(msg.contains("gnomad_genes"));
    assert!(msg.contains("clinvar_protein"));
}

#[test]
fn annotate_vcf_emits_spliceai_from_fastsa() {
    let tmp = tempfile::tempdir().unwrap();
    let spliceai_source = tmp.path().join("spliceai-mini.vcf");
    let input_vcf = tmp.path().join("input-no-spliceai-info.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("spliceai-mini");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&spliceai_source, SPLICEAI_SOURCE_VCF).unwrap();
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    run_sa_build(
        "spliceai",
        spliceai_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: Some(gff3.to_string_lossy().into_owned()),
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();

    assert!(
        annotated.contains("##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAI annotations. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">"),
        "VCF output should declare the SpliceAI INFO header:\n{}",
        annotated
    );
    assert!(
        annotated.contains("SpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13"),
        "VCF output should emit SpliceAI from fastSA for single-alt records:\n{}",
        annotated
    );
    assert!(
        annotated.contains(
            "SpliceAI=T|GENE2|0.00|0.10|0.00|0.92|3|-5|10|-2,A|GENE2|0.50|0.00|0.00|0.00|7|0|0|0"
        ),
        "VCF output should emit one SpliceAI value per matching alternate allele:\n{}",
        annotated
    );
}

#[test]
fn annotate_vcf_emits_spliceai_for_uploaded_indel_alleles() {
    let tmp = tempfile::tempdir().unwrap();
    let spliceai_source = tmp.path().join("spliceai-indel-mini.vcf");
    let input_vcf = tmp.path().join("input-no-spliceai-info.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("spliceai-indel-mini");
    let phylop_source = tmp.path().join("phylop-indel.tsv");
    let phylop_output_base = tmp.path().join("phylop-indel");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&spliceai_source, SPLICEAI_INDEL_SOURCE_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();
    fs::write(
        &phylop_source,
        "\
chr1\t26001\t3.14
chr1\t26011\t2.71
",
    )
    .unwrap();
    fs::write(
        &input_vcf,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t26000\t.\tGA\tG\t.\t.\t.
1\t26010\t.\tG\tGA\t.\t.\t.
",
    )
    .unwrap();

    run_sa_build(
        "spliceai",
        spliceai_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();
    run_sa_build(
        "phylop",
        phylop_source.to_str().unwrap(),
        phylop_output_base.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: Some(gff3.to_string_lossy().into_owned()),
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();
    assert!(
        annotated.contains("SpliceAI=G|GENE1|0.10|0.00|0.00|0.00|4|7|27|17"),
        "deletion should use uploaded ALT allele from SpliceAI source:\n{}",
        annotated
    );
    assert!(
        annotated.contains("SpliceAI=GA|GENE1|0.00|0.20|0.00|0.00|1|7|27|17"),
        "insertion should use uploaded ALT allele from SpliceAI source:\n{}",
        annotated
    );
    assert!(
        annotated.contains("FV_PHYLOP=G|3.14"),
        "positional scores should still query fastVEP's normalized indel position:\n{}",
        annotated
    );
    assert!(!annotated.contains("SpliceAI=-|"), "{annotated}");
}

#[test]
fn annotate_vcf_replaces_existing_fastvep_info() {
    let tmp = tempfile::tempdir().unwrap();
    let spliceai_source = tmp.path().join("spliceai-mini.vcf");
    let input_vcf = tmp.path().join("input-existing-info.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("spliceai-mini");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&spliceai_source, SPLICEAI_SOURCE_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();
    fs::write(
        &input_vcf,
        "\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"stale CSQ\">
##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"stale SpliceAI\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\t.\tA\tG\t.\t.\tDP=12;CSQ=old;SpliceAI=old
",
    )
    .unwrap();

    run_sa_build(
        "spliceai",
        spliceai_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: Some(gff3.to_string_lossy().into_owned()),
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();
    assert_eq!(annotated.matches("##INFO=<ID=CSQ,").count(), 1, "{annotated}");
    assert_eq!(
        annotated.matches("##INFO=<ID=SpliceAI,").count(),
        1,
        "{annotated}"
    );
    assert!(annotated.contains("DP=12;SpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13;CSQ=G|"));
    assert!(!annotated.contains("CSQ=old"), "{annotated}");
    assert!(!annotated.contains("SpliceAI=old"), "{annotated}");
    assert!(!annotated.contains("stale CSQ"), "{annotated}");
    assert!(!annotated.contains("stale SpliceAI"), "{annotated}");
}

#[test]
fn annotate_vcf_emits_fastsa_projection_for_gnomad() {
    let tmp = tempfile::tempdir().unwrap();
    let gnomad_source = tmp.path().join("gnomad-mini.vcf");
    let input_vcf = tmp.path().join("input.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("gnomad-mini");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&gnomad_source, GNOMAD_SOURCE_VCF).unwrap();
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    run_sa_build(
        "gnomad",
        gnomad_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: Some(gff3.to_string_lossy().into_owned()),
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();
    assert!(
        annotated.contains("##INFO=<ID=FV_GNOMAD,Number=.,Type=String"),
        "{annotated}"
    );
    assert!(
        annotated.contains("FV_GNOMAD=G|0.00012|12|100000|0|0.00021"),
        "{annotated}"
    );
    assert!(!annotated.contains("FV_GNOMAD={"), "{annotated}");
}
