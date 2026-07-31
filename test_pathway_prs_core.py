from pathlib import Path

import pandas as pd
import pathway_pgs_app
import ipywidgets as widgets

from pathway_prs_core import (
    WorkflowConfig,
    build_snp_set_from_variant_gene_mapping,
    build_pathway_variant_union_bed,
    build_variant_mapping_union_bed,
    build_prset_command,
    command_as_shell,
    generate_pgen_adapter_script,
    harmonize_gwas_ids_to_bed,
    harmonize_gwas_ids_to_pgen,
    infer_gwas_schema,
    plot_pathway_pgs_results,
    read_gmt_summary,
    summarize_aggregate_results,
)


def test_gmt_and_command(tmp_path: Path):
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text(
        "P1\tplain text description\tG1\tG2\n"
        "P2\thttps://example.org\tG2\tG3\tG4\n"
    )
    summary, preview = read_gmt_summary(gmt)
    assert summary["n_genes"].tolist() == [2, 3]
    assert len(preview) == 5

    config = WorkflowConfig(
        base_gwas="base.tsv",
        target_prefix="target",
        gtf_file="genes.gtf",
        gmt_file=str(gmt),
        prsice_r="PRSice.R",
        prsice_binary="PRSice_linux",
        pvalue_thresholds=[0.05, 1.0],
    )
    command = build_prset_command(config)
    assert "--msigdb" in command
    assert "--bar-levels" in command
    assert "--fastscore" in command
    assert "0.05,1" in command
    assert "--no-regress" in command
    assert "--all-score" in command
    assert "--keep" not in command
    assert "PRSice_linux" in command_as_shell(command)

    no_full = build_prset_command(
        WorkflowConfig(
            base_gwas="base.tsv",
            target_prefix="target",
            gtf_file="genes.gtf",
            gmt_file=str(gmt),
            prsice_r="PRSice.R",
            prsice_binary="PRSice_linux",
            pvalue_thresholds=[0.05],
        )
    )
    assert "--no-full" in no_full


def test_command_includes_target_keep_file():
    command = build_prset_command(
        WorkflowConfig(
            base_gwas="base.tsv",
            target_prefix="target_chr#",
            target_keep_file="participants.keep",
            gtf_file="genes.gtf",
            gmt_file="pathways.gmt",
            prsice_r="PRSice.R",
            prsice_binary="PRSice_linux",
        )
    )
    keep_index = command.index("--keep")
    assert command[keep_index + 1] == "participants.keep"


def test_pgen_adapter_uses_union_extract():
    script = generate_pgen_adapter_script(
        pgen_pattern="/data/chr{chr}",
        chromosomes=[1, 2],
        extract_file="union.txt",
        keep_file="people.keep",
        output_dir="converted",
        extract_format="bed1",
        set_coordinate_allele_ids=True,
    )
    assert script.count("--pfile") == 2
    assert script.count("--extract bed1 union.txt") == 2
    assert "--keep people.keep" in script
    assert "--set-all-var-ids '@:#:$r:$a'" in script
    assert "target_prefixes.txt" in script


def test_infer_finngen_gwas_schema(tmp_path: Path):
    gwas = tmp_path / "finngen.tsv"
    gwas.write_text(
        "#chrom\tpos\tref\talt\trsids\tpval\tbeta\taf_alt\n"
        "1\t100\tG\tA\trs1\t0.01\t0.2\t0.25\n"
    )
    schema = infer_gwas_schema(gwas)
    assert schema["separator"] == "tab"
    assert schema["statistic_type"] == "BETA"
    assert schema["gwas_columns"] == {
        "chr": "#chrom",
        "bp": "pos",
        "snp": "rsids",
        "a1": "alt",
        "a2": "ref",
        "p": "pval",
        "stat": "beta",
    }
    assert schema["maf_column"] == ""
    assert schema["info_column"] == ""
    assert schema["missing_required_fields"] == []


def test_compact_app_launches_without_writing(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(pathway_pgs_app, "HOME", tmp_path)
    monkeypatch.setattr(pathway_pgs_app, "ANALYSIS", tmp_path / "analysis")
    monkeypatch.setattr(pathway_pgs_app, "DATA", tmp_path / "analysis" / "data")
    monkeypatch.setattr(
        pathway_pgs_app,
        "RESULTS",
        tmp_path / "analysis" / "results" / "pathway_prs_tutorial",
    )
    monkeypatch.setattr(
        pathway_pgs_app,
        "WORK",
        tmp_path / "analysis" / "data" / "pathway_pgs_workspace",
    )
    app = pathway_pgs_app.launch_pathway_pgs_app()
    assert app.__class__.__name__ == "VBox"
    assert app.layout.display != "none"
    assert app.layout.width == "100%"
    assert not (tmp_path / "analysis").exists()

    pending = [app]
    observed = []
    while pending:
        widget = pending.pop()
        observed.append(widget)
        pending.extend(getattr(widget, "children", ()))
    assert not any(isinstance(widget, widgets.FileUpload) for widget in observed)
    assert any(
        "Scoring engine not found" in getattr(widget, "value", "")
        for widget in observed
        if isinstance(widget, widgets.HTML)
    )

    run_mode = next(
        widget
        for widget in observed
        if isinstance(widget, widgets.ToggleButtons)
        and widget.description == "Run mode"
    )
    run_mode.value = "full"
    assert app.layout.display != "none"


def test_build_pathway_variant_union_bed(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        '1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG1.2"; gene_name "G1";\n'
        '1\ttest\tgene\t500\t600\t.\t-\t.\tgene_id "ENSG2"; gene_name "G2";\n'
    )
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("P1\tdescription\tG1\n")
    gwas = tmp_path / "gwas.tsv"
    pd.DataFrame(
        {
            "CHR": [1, 1, 1, 1],
            "BP": [79, 80, 205, 211],
            "SNP": ["a", "b", "c", "d"],
            "A1": ["A"] * 4,
            "A2": ["G"] * 4,
            "BETA": [0.1] * 4,
            "P": [0.5] * 4,
        }
    ).to_csv(gwas, sep="\t", index=False)
    config = WorkflowConfig(
        base_gwas=str(gwas),
        gtf_file=str(gtf),
        gmt_file=str(gmt),
        window_5_bp=20,
        window_3_bp=10,
    )
    output = tmp_path / "union.bed1"
    audit = build_pathway_variant_union_bed(config, output, chunksize=2)
    observed = pd.read_csv(
        output, sep="\t", header=None, names=["chromosome", "start", "end"]
    )
    assert observed["start"].tolist() == [80, 205]
    assert audit["unique_gwas_positions_in_pathway_windows"] == 2


def test_harmonize_gwas_ids_to_pgen(tmp_path: Path):
    prefix = tmp_path / "target_chr1"
    (tmp_path / "target_chr1.pvar").write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "chr1\t80\told1\tG\tA\n"
        "chr1\t205\told2\tC\tT\n"
    )
    union = tmp_path / "union.bed1"
    union.write_text("1\t80\t80\n1\t205\t205\n")
    gwas = tmp_path / "gwas.tsv"
    pd.DataFrame(
        {
            "#CHROM": [1, 1, 1],
            "BP": [80, 205, 300],
            "SNP": ["rs1", "rs2", "rs3"],
            "A1": ["A", "T", "A"],
            "A2": ["G", "C", "C"],
            "BETA": [0.1, 0.2, 0.3],
            "P": [0.01, 0.02, 0.03],
        }
    ).to_csv(gwas, sep="\t", index=False)
    config = WorkflowConfig(
        base_gwas=str(gwas),
        gwas_columns={
            "chr": "#CHROM",
            "bp": "BP",
            "snp": "SNP",
            "a1": "A1",
            "a2": "A2",
            "stat": "BETA",
            "p": "P",
        },
    )
    output = tmp_path / "harmonized.tsv.gz"
    audit = harmonize_gwas_ids_to_pgen(
        config,
        pgen_pattern=str(tmp_path / "target_chr{chr}"),
        chromosomes=[1],
        candidate_bed1_file=union,
        output_gwas=output,
        chunksize=2,
    )
    observed = pd.read_csv(output, sep="\t")
    assert observed.columns.tolist() == [
        "CHR", "BP", "SNP", "A1", "A2", "BETA", "P", "ORIGINAL_SNP_ID"
    ]
    assert observed["SNP"].tolist() == ["chr1:80:G:A", "chr1:205:C:T"]
    assert observed["ORIGINAL_SNP_ID"].tolist() == ["rs1", "rs2"]
    assert audit["unique_harmonized_target_ids"] == 2
    safe_schema = infer_gwas_schema(output)
    safe_command = build_prset_command(
        WorkflowConfig(
            base_gwas=str(output),
            target_prefix="target",
            gtf_file="genes.gtf",
            gmt_file="pathways.gmt",
            prsice_r="PRSice.R",
            prsice_binary="PRSice_linux",
            gwas_columns=safe_schema["gwas_columns"],
        )
    )
    assert safe_command[safe_command.index("--chr") + 1] == "CHR"


def test_harmonize_gwas_ids_to_existing_bed(tmp_path: Path):
    prefix = tmp_path / "target_chr1"
    (tmp_path / "target_chr1.bim").write_text(
        "1\taou_variant_1\t0\t80\tA\tG\n"
        "1\taou_variant_2\t0\t205\tT\tC\n"
    )
    for suffix in (".bed", ".fam"):
        (tmp_path / f"target_chr1{suffix}").write_text("")
    union = tmp_path / "union.bed1"
    union.write_text("1\t80\t80\n1\t205\t205\n")
    gwas = tmp_path / "gwas.tsv"
    pd.DataFrame(
        {
            "CHR": [1, 1],
            "BP": [80, 205],
            "SNP": ["rs1", "rs2"],
            "A1": ["A", "T"],
            "A2": ["G", "C"],
            "BETA": [0.1, 0.2],
            "P": [0.01, 0.02],
        }
    ).to_csv(gwas, sep="\t", index=False)
    config = WorkflowConfig(
        base_gwas=str(gwas),
        target_prefix=str(tmp_path / "target_chr#"),
    )
    output = tmp_path / "harmonized.tsv.gz"
    audit = harmonize_gwas_ids_to_bed(
        config,
        candidate_bed1_file=union,
        output_gwas=output,
        chunksize=1,
    )
    observed = pd.read_csv(output, sep="\t")
    assert observed.columns.tolist() == [
        "CHR", "BP", "SNP", "A1", "A2", "BETA", "P", "ORIGINAL_SNP_ID"
    ]
    assert observed["SNP"].tolist() == ["aou_variant_1", "aou_variant_2"]
    assert observed["ORIGINAL_SNP_ID"].tolist() == ["rs1", "rs2"]
    assert audit["unique_harmonized_target_ids"] == 2


def test_command_rejects_hash_prefixed_gwas_header():
    config = WorkflowConfig(
        base_gwas="base.tsv",
        target_prefix="target",
        gtf_file="genes.gtf",
        gmt_file="pathways.gmt",
        prsice_r="PRSice.R",
        prsice_binary="PRSice_linux",
        gwas_columns={
            "chr": "#chrom",
            "bp": "pos",
            "snp": "rsids",
            "a1": "alt",
            "a2": "ref",
            "stat": "beta",
            "p": "pval",
        },
    )
    try:
        build_prset_command(config)
    except ValueError as error:
        assert "Run Step 2" in str(error)
    else:
        raise AssertionError("Unsafe hash-prefixed header was accepted")


def test_custom_variant_gene_mapping_builds_pathway_snp_sets(tmp_path: Path):
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "rsid\tgene\n"
        "rs1\tENSG1.2\n"
        "rs2\tENSG1\n"
        "rs3\tENSG2\n"
    )
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("P1\tdescription\tENSG1\nP2\tdescription\tENSG2\n")
    harmonized = tmp_path / "harmonized.tsv"
    harmonized.write_text(
        "SNP\tORIGINAL_SNP_ID\n"
        "aou1\trs1\n"
        "aou2\trs2\n"
        "aou3\trs3\n"
    )
    output = tmp_path / "sets.gmt"
    audit = build_snp_set_from_variant_gene_mapping(
        mapping_file=mapping,
        gmt_file=gmt,
        output_file=output,
        harmonized_gwas_file=harmonized,
        harmonized_snp_column="SNP",
    )
    lines = output.read_text().splitlines()
    assert lines[0].split("\t")[2:] == ["aou1", "aou2"]
    assert lines[1].split("\t")[2:] == ["aou3"]
    assert audit["pathways_with_variants"] == 2


def test_coordinate_gene_mapping_builds_union_and_pathway_sets(tmp_path: Path):
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "chromosome\tposition\tgene\n"
        "1\t80\tENSG1\n"
        "1\t205\tENSG2\n"
    )
    gwas = tmp_path / "gwas.tsv"
    pd.DataFrame(
        {
            "CHR": [1, 1],
            "BP": [80, 205],
            "SNP": ["rs1", "rs2"],
            "A1": ["A", "T"],
            "A2": ["G", "C"],
            "BETA": [0.1, 0.2],
            "P": [0.01, 0.02],
        }
    ).to_csv(gwas, sep="\t", index=False)
    union = tmp_path / "union.bed1"
    audit = build_variant_mapping_union_bed(
        WorkflowConfig(base_gwas=str(gwas)), mapping, union
    )
    assert audit["mapping_mode"] == "chromosome_position"
    assert union.read_text().splitlines() == ["1\t80\t80", "1\t205\t205"]

    harmonized = tmp_path / "harmonized.tsv"
    pd.DataFrame(
        {
            "CHR": [1, 1],
            "BP": [80, 205],
            "SNP": ["aou1", "aou2"],
            "ORIGINAL_SNP_ID": ["rs1", "rs2"],
        }
    ).to_csv(harmonized, sep="\t", index=False)
    gmt = tmp_path / "pathways.gmt"
    gmt.write_text("P1\tdescription\tENSG1\nP2\tdescription\tENSG2\n")
    output = tmp_path / "sets.gmt"
    set_audit = build_snp_set_from_variant_gene_mapping(
        mapping_file=mapping,
        gmt_file=gmt,
        output_file=output,
        harmonized_gwas_file=harmonized,
        harmonized_snp_column="SNP",
    )
    assert set_audit["pathways_with_variants"] == 2
    assert output.read_text().splitlines() == [
        "P1\tdescription\taou1",
        "P2\tdescription\taou2",
    ]


def test_score_only_results_are_aggregated_without_ids(tmp_path: Path):
    score_file = tmp_path / "demo.all_score"
    score_file.write_text(
        "FID IID Pathway_A Pathway_B\n"
        "1 10 0.1 -0.2\n"
        "2 20 0.3 0.4\n"
    )
    summary = summarize_aggregate_results(
        WorkflowConfig(output_dir=str(tmp_path), output_prefix="demo")
    )
    assert summary["score_name"].tolist() == ["Pathway_A", "Pathway_B"]
    assert summary["n"].tolist() == [2, 2]
    assert "IID" not in summary.columns
    assert (tmp_path / "demo.aggregate_score_summary.tsv").exists()
    figure = plot_pathway_pgs_results(summary)
    assert len(figure.axes) == 2


def test_example_config_is_loadable():
    source = Path(__file__).resolve().parent / "example_config.json"
    config = WorkflowConfig.from_json(source)
    assert config.genome_build == "GRCh38"
    assert config.pathway_input_mode == "gtf_gmt"
