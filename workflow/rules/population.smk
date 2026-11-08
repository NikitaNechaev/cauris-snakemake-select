rule genome_qc:
    input:
        code=script("cauris_population_stage"),
        selected=SRC["selected"],
        genomes=genome_inputs,
    output:
        directory(result(GENOME_QC)),
    params:
        genomes=SRC["genomes"],
    log:
        workflow_log("genome_qc"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} qc --selected {input.selected:q} --genomes {params.genomes:q} --qc {output:q} "
        "> {log:q} 2>&1"

if not SRC["ani_matrix"]:
    rule skani:
        input:
            qc=result(GENOME_QC),
            genomes=genome_inputs,
        output:
            ani=ensure(ANI_MATRIX, non_empty=True),
            af=ensure(AF_MATRIX, non_empty=True),
        threads: T["search"]
        log:
            workflow_log("skani"),
        conda:
            "../envs/skani.yaml"
        shell:
            "skani triangle --full-matrix -l {input.qc:q}/genome_paths.txt "
            "-t {threads} -o {output.ani:q} > {log:q} 2>&1"

rule ani_qc:
    input:
        code=script("analyze_cauris_ani"),
        qc=result(GENOME_QC),
        ani=ANI_MATRIX,
        af=AF_MATRIX,
    output:
        directory(result(ANI)),
    params:
        clone=A["clone_ani"],
        group=A["group_ani"],
        af=A["min_af"],
    log:
        workflow_log("ani_qc"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --ani {input.ani:q} --af {input.af:q} --manifest {input.qc:q}/genome_manifest.tsv "
        "--clone-ani {params.clone} --group-ani {params.group} --min-af {params.af} --out {output:q} > {log:q} 2>&1"

rule clade_tables:
    input:
        code=script("prepare_inputs"),
        population=script("cauris_population_stage"),
        qc=result(GENOME_QC),
        assignments=SRC["clade_assignments"],
        comparisons=SRC["clade_comparisons"],
    output:
        assignments=result(f"{CLADES}/clade_assignments.tsv"),
        comparisons=result(f"{CLADES}/all_reference_comparisons.tsv"),
        summary=result(f"{CLADES}/clade_summary.json"),
    log:
        workflow_log("clade_tables"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} clades --manifest {input.qc:q}/genome_manifest.tsv "
        "--assignments {input.assignments:q} --comparisons {input.comparisons:q} "
        "--summary {output.summary:q} > {log:q} 2>&1"

rule clade_diagnostics:
    input:
        code=script("diagnose_cauris_clades"),
        assignments=result(f"{CLADES}/clade_assignments.tsv"),
        comparisons=result(f"{CLADES}/all_reference_comparisons.tsv"),
    output:
        result(f"{CLADES}/clade_diagnostics.tsv"),
    params:
        folder=lambda wc, output: str(Path(output[0]).parent),
        af=A["min_af"],
    log:
        workflow_log("clade_diagnostics"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --folder {params.folder:q} --min-af {params.af} > {log:q} 2>&1"

rule map_genome:
    input:
        code=script("cauris_population_stage"),
        qc=result(GENOME_QC),
        selected=SRC["selected"],
        genome=lambda wc: assembly_genome(wc.assembly),
        candidates=result(f"{CANDIDATES}/auris_priority.faa"),
    output:
        gff=ensure(result(f"{MAPPING}/{{assembly}}.gff"), non_empty=True),
        stamp=ensure(result(f"{MAPPING}/{{assembly}}.done.json"), non_empty=True),
    params:
        folder=lambda wc, output: str(Path(output.gff).parent),
    threads: T["mapping"]
    log:
        run=workflow_log("mapping/{assembly}"),
        aligner=result(f"{MAPPING}/{{assembly}}.log"),
    conda:
        "../envs/miniprot.yaml"
    shell:
        "python {input.code:q} map --manifest {input.qc:q}/genome_manifest.tsv --selected {input.selected:q} "
        "--candidates {input.candidates:q} --alignments {params.folder:q} --assembly {wildcards.assembly:q} "
        "--code 12 --threads {threads} --force > {log.run:q} 2>&1"

rule mapping_summary:
    input:
        code=script("cauris_population_stage"),
        qc=result(GENOME_QC),
        selected=SRC["selected"],
        candidates=result(f"{CANDIDATES}/auris_priority.faa"),
        gff=lambda wc: per_assembly(wc, f"{MAPPING}/{{assembly}}.gff"),
        stamps=lambda wc: per_assembly(wc, f"{MAPPING}/{{assembly}}.done.json"),
    output:
        details=result(f"{MAP_SUMMARY}/candidate_genome_details.tsv"),
        binary=result(f"{MAP_SUMMARY}/candidate_presence_matrix.tsv"),
        status=result(f"{MAP_SUMMARY}/candidate_alignment_status_matrix.tsv"),
        counts=result(f"{MAP_SUMMARY}/candidate_presence_summary.tsv"),
        meta=result(f"{MAP_SUMMARY}/summary.json"),
    params:
        alignments=result(MAPPING),
        out=lambda wc, output: str(Path(output.details).parent),
        identity=A["mapping_identity"],
        coverage=A["mapping_coverage"],
    log:
        workflow_log("mapping_summary"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} summarize --manifest {input.qc:q}/genome_manifest.tsv --selected {input.selected:q} "
        "--candidates {input.candidates:q} --alignments {params.alignments:q} --results {params.out:q} "
        "--min-identity {params.identity} --min-coverage {params.coverage} > {log:q} 2>&1"

rule mapping_audit:
    input:
        code=script("audit_miniprot_by_clade"),
        qc=result(GENOME_QC),
        genomes=genome_inputs,
        details=result(f"{MAP_SUMMARY}/candidate_genome_details.tsv"),
        meta=result(f"{MAP_SUMMARY}/summary.json"),
        clades=result(f"{CLADES}/clade_assignments.tsv"),
        annotations=result(f"{HUMAN}/auris_cluster_audit.tsv"),
    output:
        directory(result(AUDIT)),
    params:
        window=A["edge_window"],
        min_n=A["min_clade_n"],
        identity=A["near_identity"],
        coverage=A["near_coverage"],
    log:
        workflow_log("mapping_audit"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --manifest {input.qc:q}/genome_manifest.tsv --details {input.details:q} "
        "--clades {input.clades:q} --metadata {input.annotations:q} --out {output:q} "
        "--edge-check --edge-window {params.window} --min-clade-n {params.min_n} "
        "--near-identity {params.identity} --near-coverage {params.coverage} > {log:q} 2>&1"

rule legacy_genome_qc:
    input:
        code=script("cauris_population_stage_legacy"),
        selected=SRC["selected"],
        genomes=genome_inputs,
    output:
        directory(result(f"{LEGACY}/qc")),
    params:
        genomes=SRC["genomes"],
    log:
        workflow_log("legacy_genome_qc"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} qc --selected {input.selected:q} --genomes {params.genomes:q} "
        "--out {output:q} > {log:q} 2>&1"

rule legacy_mapping_summary:
    input:
        code=script("cauris_population_stage_legacy"),
        qc=result(f"{LEGACY}/qc"),
        gff=lambda wc: per_assembly(wc, f"{MAPPING}/{{assembly}}.gff"),
        stamps=lambda wc: per_assembly(wc, f"{MAPPING}/{{assembly}}.done.json"),
        candidates=result(f"{CANDIDATES}/auris_priority.faa"),
    output:
        directory(result(f"{LEGACY}/summary")),
    params:
        alignments=result(MAPPING),
        identity=A["mapping_identity"],
        coverage=A["mapping_coverage"],
    log:
        workflow_log("legacy_mapping_summary"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} summarize --manifest {input.qc:q}/genome_manifest.tsv "
        "--gff-dir {params.alignments:q} --stamp-suffix .done.json --candidates {input.candidates:q} "
        "--min-identity {params.identity} --min-coverage {params.coverage} --out {output:q} > {log:q} 2>&1"
