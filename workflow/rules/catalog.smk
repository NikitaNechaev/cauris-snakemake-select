if CAT["mode"] == "snapshot":
    rule catalog:
        input:
            code=script("prepare_inputs"),
            snapshots=[str(Path(CAT["directory"]) / f) for f in CATALOG_FILES],
            optional=[str(Path(CAT["directory"]) / f) for f in CATALOG_OPTIONAL
                      if (Path(CAT["directory"]) / f).is_file()],
        output:
            directory(result(CATALOG)),
        log:
            workflow_log("catalog"),
        conda:
            "../envs/python.yaml"
        shell:
            "python {input.code:q} snapshot --out {output:q} {input.snapshots:q} {input.optional:q} > {log:q} 2>&1"
else:
    rule catalog:
        input:
            code=script("download_surface_antigens_v2"),
        output:
            directory(result(CATALOG)),
        params:
            scope=CAT["scope"],
            delay=CAT["request_delay"],
        log:
            workflow_log("catalog"),
        conda:
            "../envs/python.yaml"
        shell:
            "python {input.code:q} --scope {params.scope:q} --sleep {params.delay} --out {output:q} > {log:q} 2>&1"

rule catalog_qc:
    input:
        code=script("qc_surface_catalog"),
        catalog=result(CATALOG),
    output:
        directory(result(QC)),
    log:
        workflow_log("catalog_qc"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --input {input.catalog:q} --out {output:q} > {log:q} 2>&1"

rule cluster_catalog:
    input:
        qc=result(QC),
    output:
        result(f"{CLUSTERS}/surface_cluster.tsv"),
    params:
        identity=A["cluster_identity"],
        coverage=A["cluster_coverage"],
    threads: T["search"]
    log:
        workflow_log("mmseqs")
    conda:
        "../envs/mmseqs.yaml"
    shell:
        r"""
        tmp=$(mktemp -d)
        trap 'rm -rf -- "$tmp"' EXIT
        mmseqs easy-cluster {input.qc:q}/surface_priority.faa "$tmp/surface" "$tmp/work" \
            --min-seq-id {params.identity} -c {params.coverage} --cov-mode 0 \
            --threads {threads} > {log:q} 2>&1
        mv "$tmp/surface_cluster.tsv" {output:q}
        """

rule cluster_summary:
    input:
        code=script("summarize_surface_clusters"),
        clusters=result(f"{CLUSTERS}/surface_cluster.tsv"),
        qc=result(QC),
    output:
        result(f"{CLUSTERS}/cluster_summary.tsv"),
    log:
        workflow_log("cluster_summary"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --clusters {input.clusters:q} --metadata {input.qc:q}/surface_priority.tsv "
        "--out {output:q} > {log:q} 2>&1"

rule select_candidates:
    input:
        code=script("prepare_inputs"),
        qc=result(QC),
        clusters=result(f"{CLUSTERS}/surface_cluster.tsv"),
        summary=result(f"{CLUSTERS}/cluster_summary.tsv"),
        allowlist=[SRC["candidate_accessions"]] if SRC["candidate_accessions"] else [],
    output:
        faa=ensure(result(f"{CANDIDATES}/auris_priority.faa"), non_empty=True),
        metadata=ensure(result(f"{HUMAN}/auris_cluster_audit.tsv"), non_empty=True),
    params:
        allowlist_flag="--allowlist" if SRC["candidate_accessions"] else "",
    log:
        workflow_log("select_candidates"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} candidates --qc {input.qc:q} --clusters {input.clusters:q} "
        "--fasta {output.faa:q} --metadata {output.metadata:q} "
        "{params.allowlist_flag} {input.allowlist:q} > {log:q} 2>&1"

rule human_database:
    input:
        SRC["human_proteome"],
    output:
        result(f"{HUMAN}/human.dmnd"),
    log:
        workflow_log("human_database"),
    conda:
        "../envs/diamond.yaml"
    shell:
        "diamond makedb --in {input:q} --db {output:q} --threads {threads} > {log:q} 2>&1"

rule human_search:
    input:
        query=result(f"{CANDIDATES}/auris_priority.faa"),
        db=result(f"{HUMAN}/human.dmnd"),
    output:
        result(f"{HUMAN}/human_hits.tsv"),
    threads: T["search"]
    params:
        fields=HIT_FIELDS,
        evalue=A["diamond_evalue"],
        targets=A["diamond_max_targets"],
    log:
        workflow_log("human_search"),
    conda:
        "../envs/diamond.yaml"
    shell:
        "diamond blastp --query {input.query:q} --db {input.db:q} --out {output:q} --threads {threads} "
        "--sensitive --evalue {params.evalue} --max-target-seqs {params.targets} "
        "--outfmt 6 {params.fields} > {log:q} 2>&1"

rule human_summary:
    input:
        code=script("summarize_human_similarity"),
        hits=result(f"{HUMAN}/human_hits.tsv"),
        metadata=result(f"{HUMAN}/auris_cluster_audit.tsv"),
    output:
        result(f"{HUMAN}/auris_human_summary.tsv"),
    log:
        workflow_log("human_summary"),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.code:q} --hits {input.hits:q} --metadata {input.metadata:q} --out {output:q} > {log:q} 2>&1"

if SRC["isolates"]:
    rule isolate_prepare:
        input:
            code=script("check_auris_isolates"),
            proteomes=isolate_inputs,
        output:
            directory(result(f"{ISOLATES}/prepared")),
        params:
            directory=SRC["isolates"],
        log:
            workflow_log("isolate_prepare"),
        conda:
            "../envs/python.yaml"
        shell:
            "python {input.code:q} prepare --ncbi-dir {params.directory:q} --out {output:q} > {log:q} 2>&1"

    rule isolate_database:
        input:
            result(f"{CANDIDATES}/auris_priority.faa"),
        output:
            result(f"{ISOLATES}/candidates.dmnd"),
        log:
            workflow_log("isolate_database"),
        conda:
            "../envs/diamond.yaml"
        shell:
            "diamond makedb --in {input:q} --db {output:q} --threads {threads} > {log:q} 2>&1"

    rule isolate_search:
        input:
            prepared=result(f"{ISOLATES}/prepared"),
            db=result(f"{ISOLATES}/candidates.dmnd"),
        output:
            result(f"{ISOLATES}/hits.tsv"),
        params:
            evalue=A["diamond_evalue"],
            fields=HIT_FIELDS,
        threads: T["search"]
        log:
            workflow_log("isolate_search"),
        conda:
            "../envs/diamond.yaml"
        shell:
            "diamond blastp --query {input.prepared:q}/combined_proteomes.faa --db {input.db:q} "
            "--out {output:q} --threads {threads} --sensitive --evalue {params.evalue} "
            "--max-target-seqs 0 --outfmt 6 {params.fields} > {log:q} 2>&1"

    rule isolate_summary:
        input:
            code=script("check_auris_isolates"),
            hits=result(f"{ISOLATES}/hits.tsv"),
            candidates=result(f"{CANDIDATES}/auris_priority.faa"),
            prepared=result(f"{ISOLATES}/prepared"),
        output:
            directory(result(f"{ISOLATES}/summary")),
        params:
            identity=A["isolate_identity"],
            protein_coverage=A["isolate_protein_coverage"],
            candidate_coverage=A["isolate_candidate_coverage"],
        log:
            workflow_log("isolate_summary"),
        conda:
            "../envs/python.yaml"
        shell:
            "python {input.code:q} summarize --hits {input.hits:q} --candidates {input.candidates:q} "
            "--counts {input.prepared:q}/assembly_counts.tsv --out {output:q} --min-id {params.identity} "
            "--min-protein-cov {params.protein_coverage} --min-candidate-cov {params.candidate_coverage} > {log:q} 2>&1"
