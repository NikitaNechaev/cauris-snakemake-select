import json
import re
import sys
from pathlib import Path

from snakemake.utils import min_version, validate

min_version("9.23.0")
configfile: "config/config.yaml"
validate(config, str(WORKFLOW_DIR.parent / "config/schema.yaml"))
SCRIPTS = WORKFLOW_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS))
from cauris_population_stage import FASTA_SUFFIXES, read_selected

ROOT = Path(config["results"]).resolve()
SRC = config["inputs"]
CAT = config["catalog"]
A = config["analysis"]
T = config["threads"]

def result(path):
    return str(ROOT / path)

def script(name):
    return str(SCRIPTS / f"{name}.py")

def workflow_log(name):
    return result(f"logs/{name}.log")

def assemblies():
    return read_selected(SRC["selected"])

def assembly_genome(assembly):
    folder = Path(SRC["genomes"]) / assembly
    files = sorted(p for p in folder.glob("*") if p.is_file() and p.name.lower().endswith(FASTA_SUFFIXES))
    if len(files) != 1:
        raise ValueError(f"Expected one genomic FASTA in {folder}, found {len(files)}")
    return str(files[0].resolve())

def genome_inputs(wildcards):
    return [assembly_genome(a) for a in assemblies()]

def per_assembly(wildcards, pattern):
    return expand(result(pattern), assembly=assemblies())

def isolate_inputs(wildcards):
    files = sorted(Path(SRC["isolates"]).rglob("protein.faa"))
    files = [str(p.resolve()) for p in files if ASSEMBLY_RE.fullmatch(p.parent.name)]
    if not files:
        raise ValueError(f"No accession/protein.faa files under {SRC['isolates']}")
    return files

ASSEMBLY_RE = re.compile(r"GC[AF]_\d+\.\d+")
HIT_FIELDS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
CATALOG = "01_surface_antigens/data"
CATALOG_FILES = ("surface_candidates.faa", "surface_candidates.tsv", "run_metadata.json")
CATALOG_OPTIONAL = ("query_audit.tsv", "by_organism.tsv", "cell_wall_annotated.faa")
QC = "02_qc_catalog"
CLUSTERS = "03_clusters"
CANDIDATES = "04_summarize_clusters"
HUMAN = "06_blastp_on_human"
ISOLATES = "07_annotated_isolates"
POP = "08_population"
GENOME_QC = f"{POP}/qc"
MAPPING = f"{POP}/alignments"
MAP_SUMMARY = f"{POP}/results"
ANI = f"{MAP_SUMMARY}/ani_population"
CLADES = f"{MAP_SUMMARY}/clade_typing"
AUDIT = f"{MAP_SUMMARY}/mapping_audit"
PROTEINS = "09_protein_sequences"
MAFFT = "10_mafft"
CONSERVATION = "11_conservation_statistics"
FINAL = "12_final_report"
FINAL_FILES = (
    "integrated_candidates.tsv",
    "candidate_by_clade.tsv",
    "evidence_issues.tsv",
    "input_manifest_sha256.tsv",
    "raw_inputs_sha256.tsv",
    "raw_config.json",
    "final_summary.json",
    "RESULTS_RU.md",
    "README_RU.md",
    "figures/01_extraction_by_clade.png",
    "figures/02_recovery_vs_entropy.png",
    "figures/03_human_reported_hits.png",
    "figures/04_mapping_qc_flags.png",
    "final_results_bundle.zip",
)
LEGACY = "13_legacy_population"

if bool(SRC["ani_matrix"]) != bool(SRC["af_matrix"]):
    raise ValueError("Set both inputs.ani_matrix and inputs.af_matrix, or leave both null")
ANI_MATRIX = SRC["ani_matrix"] or result(f"{POP}/skani_matrix.txt")
AF_MATRIX = SRC["af_matrix"] or result(f"{POP}/skani_matrix.txt.af")

def provenance_inputs(wildcards):
    inputs = [
        SRC["selected"], SRC["human_proteome"], SRC["clade_assignments"], SRC["clade_comparisons"],
        *genome_inputs(wildcards),
        *[f for f in (SRC["candidate_accessions"], SRC["ani_matrix"], SRC["af_matrix"]) if f],
    ]
    if SRC["isolates"]:
        inputs.extend(isolate_inputs(wildcards))
    if CAT["mode"] == "snapshot":
        folder = Path(CAT["directory"])
        inputs.extend(str(folder / name) for name in CATALOG_FILES)
        inputs.extend(str(folder / name) for name in CATALOG_OPTIONAL if (folder / name).is_file())
    return inputs
