TRANSLATIONS = f"{POP}/translations"
PROVENANCE = "00_provenance"

rule translate_genome:
    input:
        code=script("extract_cauris_proteins"),
        qc=result(GENOME_QC),
        genome=lambda wc: assembly_genome(wc.assembly),
        gff=result(f"{MAPPING}/{{assembly}}.gff"),
        candidates=result(f"{CANDIDATES}/auris_priority.faa"),
        details=result(f"{MAP_SUMMARY}/candidate_genome_details.tsv"),
        clades=result(f"{CLADES}/clade_assignments.tsv"),
        audit=result(AUDIT),
    output:
        ensure(result(f"{TRANSLATIONS}/{{assembly}}.json"), non_empty=True),
    params:
        root=lambda wc: str(ROOT),
        cache=lambda wc, output: str(Path(output[0]).parent),
    threads: T["mapping"]
    log:
        run=workflow_log("translation/{assembly}"),
        aligner=workflow_log("translation/{assembly}.miniprot"),
    conda:
        "../envs/miniprot.yaml"
    shell:
        "python {input.code:q} --root {params.root:q} --out {params.cache:q} --cache-dir {params.cache:q} "
        "--assembly {wildcards.assembly:q} --alignment-log {log.aligner:q} --threads {threads} --force > {log.run:q} 2>&1"

rule collect_proteins:
    input:
        code=script("extract_cauris_proteins"),
        qc=result(GENOME_QC),
        translations=lambda wc: per_assembly(wc, f"{TRANSLATIONS}/{{assembly}}.json"),
        gff=lambda wc: per_assembly(wc, f"{MAPPING}/{{assembly}}.gff"),
        candidates=result(f"{CANDIDATES}/auris_priority.faa"),
        details=result(f"{MAP_SUMMARY}/candidate_genome_details.tsv"),
        clades=result(f"{CLADES}/clade_assignments.tsv"),
        audit=result(AUDIT),
    output:
        directory(result(PROTEINS)),
    params:
        root=lambda wc: str(ROOT),
        cache=lambda wc, input: str(Path(input.translations[0]).parent),
    log:
        workflow_log("collect_proteins"),
    conda:
        "../envs/miniprot.yaml"
    shell:
        "python {input.code:q} --root {params.root:q} --cache-dir {params.cache:q} --out {output:q} "
        "--collect > {log:q} 2>&1"

rule align_proteins:
    input:
        code=script("run_mafft_step19"),
        proteins=result(PROTEINS),
    output:
        directory(result(MAFFT)),
    params:
        root=lambda wc: str(ROOT),
        method=A["mafft_method"],
    threads: T["mafft"]
    log:
        workflow_log("mafft"),
    conda:
        "../envs/mafft.yaml"
    shell:
        "python {input.code:q} --root {params.root:q} --input {input.proteins:q} --out {output:q} "
        "--threads {threads} --method {params.method:q} > {log:q} 2>&1"

rule conservation:
    input:
        code=script("analyze_conservation_step20"),
        qc=result(GENOME_QC),
        clades=result(f"{CLADES}/clade_assignments.tsv"),
        proteins=result(PROTEINS),
        alignment=result(MAFFT),
    output:
        directory(result(CONSERVATION)),
    params:
        root=lambda wc: str(ROOT),
        occupancy=A["min_occupancy"],
        min_n=A["min_clade_n"],
        plot_top=A["plot_top"],
    log:
        workflow_log("conservation"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --root {params.root:q} --out {output:q} --min-occupancy {params.occupancy} "
        "--min-clade-n {params.min_n} --plot-top {params.plot_top} > {log:q} 2>&1"

rule provenance:
    input:
        code=script("prepare_inputs"),
        population=script("cauris_population_stage"),
        catalog=result(CATALOG),
        raw=provenance_inputs,
        codebase=sorted(str(p) for p in WORKFLOW_DIR.rglob("*")
                        if p.is_file() and (p.name == "Snakefile" or p.suffix in (".smk", ".py", ".yaml"))),
    output:
        manifest=result(f"{PROVENANCE}/inputs_sha256.tsv"),
        config=result(f"{PROVENANCE}/config.json"),
    params:
        configuration=json.dumps(config, sort_keys=True),
        exports=" ".join(CATALOG_FILES),
    log:
        workflow_log("provenance"),
    conda:
        "../envs/python.yaml"
    shell:
        r"""
        exports=()
        for name in {params.exports}; do exports+=({input.catalog:q}/"$name"); done
        python {input.code:q} provenance --manifest {output.manifest:q} --config-out {output.config:q} \
          --config-json {params.configuration:q} {input.raw:q} {input.codebase:q} "${{exports[@]}}" > {log:q} 2>&1
        """

rule report:
    input:
        code=script("build_final_report_step21"),
        catalog=result(CATALOG),
        catalog_qc=result(QC),
        qc=result(GENOME_QC),
        annotations=result(f"{HUMAN}/auris_cluster_audit.tsv"),
        human=result(f"{HUMAN}/auris_human_summary.tsv"),
        audit=result(AUDIT),
        ani=result(ANI),
        clades=result(f"{CLADES}/clade_assignments.tsv"),
        clade_metadata=result(f"{CLADES}/clade_summary.json"),
        diagnostics=result(f"{CLADES}/clade_diagnostics.tsv"),
        proteins=result(PROTEINS),
        alignment=result(MAFFT),
        conservation=result(CONSERVATION),
        provenance=result(f"{PROVENANCE}/inputs_sha256.tsv"),
        configuration=result(f"{PROVENANCE}/config.json"),
    output:
        files=[ensure(result(f"{FINAL}/{name}"), non_empty=True) for name in FINAL_FILES],
    params:
        root=lambda wc: str(ROOT),
        min_n=A["min_clade_n"],
        out=lambda wc, output: str(Path(output.files[0]).parent),
    log:
        workflow_log("report"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --root {params.root:q} --out {params.out:q} --min-clade-n {params.min_n} > {log:q} 2>&1"
