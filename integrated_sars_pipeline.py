#!/usr/bin/env python3
"""Integrated SARS-CoV-2 Spike Amplicon Pipeline.

Processes ARTIC amplicon sequencing data through quality filtering,
amplicon retrieval, clustering, translation, alignment, and entropy analysis.
"""

import os
import io
import re
import shutil
import subprocess
import tempfile
import argparse
import logging
from datetime import datetime
from pathlib import Path
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict, Counter
import math
from math import log

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
REF_GENOME = BASE_DIR / "WG_wuhan.fa"
SPIKE_DB = BASE_DIR / "sars_spike_database.fasta"
WUHAN_AA_SPIKE = BASE_DIR / "wuhan_AA_Spike.fasta"
AMPLICON_REF_DIR = BASE_DIR / "Amplicon_references"
PRIMER_DIR = BASE_DIR / "Primer_schemes"

MMSEQS_EXECUTABLE = "mmseqs"
MAFFT_EXECUTABLE = "mafft"
GAP_CHARS = {"-", "."}
TAIL_GAP_MIN = 2
TAIL_NON_GAP_MAX = 3

ARTIC_CONFIGS = {
    "ARTIC_3": {
        "primer_csv": PRIMER_DIR / "primers_ARTICv.3.csv",
        "positions": {
            "A71": "21450-21599", "A72": "21774-21923", "A73": "22073-22222",
            "A74": "22365-22494", "A75": "22663-22782", "A76": "22930-23079",
            "A77": "23260-23408", "A78": "23584-23733", "A79": "23890-24039",
            "A80": "24204-24353", "A81": "24505-24654", "A82": "24810-24959",
            "A83": "25103-25252", "A84": "25432-25580",
        },
        "amplicon_ref_subdir": "ARTIC3",
        "last_amplicon": "A84",
    },
    "ARTIC_4": {
        "primer_csv": PRIMER_DIR / "primers_ARTICv.4.1.csv",
        "positions": {
            "A72": "21721-21845", "A73": "21951-22079",
            "A74": "22284-22396", "A75": "22552-22700", "A76": "22814-22936",
            "A77": "23149-23213", "A78": "23376-23525", "A79": "23668-23817",
            "A80": "23981-24130", "A81": "24269-24418", "A82": "24587-24735",
            "A83": "24866-25014", "A84": "25217-25366",
        },
        "amplicon_ref_subdir": "ARTIC4",
        "last_amplicon": "A84",
    },
    "ARTIC_5": {
        "primer_csv": PRIMER_DIR / "primers_ARTICv.5.3.2.csv",
        "positions": {
            "A69": "21391-21540", "A70": "21733-21857", "A71": "21986-22135",
            "A72": "22292-22441", "A73": "22570-22719", "A74": "22902-23051",
            "A75": "23149-23223", "A76": "23489-23556", "A77": "23653-23802",
            "A78": "23980-24129", "A79": "24261-24410", "A80": "24580-24727",
            "A81": "24872-25021", "A82": "25218-25367",
        },
        "amplicon_ref_subdir": "ARTIC5",
        "last_amplicon": "A82",
    },
}

EXPECTED_FILENAME_PATTERN = re.compile(r"^[A-Za-z0-9_-]+\.fastq\.gz$")
RESERVED_SUFFIXES = ["_trimmed", "_clusters", "_curated", "_filtered", "_translated", "_aligned"]


def validate_fastq_filenames(artic_dir, artic_name):
    """Filter FASTQ files in *artic_dir* to those with valid sample names.

    A valid filename matches ``<samplename>.fastq.gz`` where *samplename*
    contains only alphanumeric characters, hyphens, and underscores and does
    not end with a reserved pipeline suffix (e.g. ``_trimmed``).

    Parameters
    ----------
    artic_dir : Path
        Directory containing ``.fastq.gz`` files for one ARTIC version.
    artic_name : str
        ARTIC version label used in log messages (e.g. ``"ARTIC_3"``).

    Returns
    -------
    list[Path]
        Sorted list of paths to FASTQ files that passed validation.
    """
    fastq_files = sorted(artic_dir.glob("*.fastq.gz"))
    valid_files = []
    for fq in fastq_files:
        fname = fq.name
        if not EXPECTED_FILENAME_PATTERN.match(fname):
            logger.warning(
                f"[{artic_name}] Skipping '{fname}': filename must match "
                f"'<samplename>.fastq.gz' where <samplename> contains only "
                f"alphanumeric characters, hyphens, and underscores. "
                f"Spaces, dots, or other special characters in the sample name "
                f"are not supported and may cause downstream steps to fail."
            )
            continue
        sample_name = fname.replace(".fastq.gz", "")
        if any(sample_name.endswith(s) for s in RESERVED_SUFFIXES):
            logger.warning(
                f"[{artic_name}] Skipping '{fname}': sample name must not end "
                f"with reserved suffixes ({', '.join(RESERVED_SUFFIXES)})."
            )
            continue
        valid_files.append(fq)
    if len(valid_files) < len(fastq_files):
        logger.warning(
            f"[{artic_name}] {len(fastq_files) - len(valid_files)} file(s) skipped "
            f"due to invalid filenames. Expected format: '<samplename>.fastq.gz' "
            f"where <samplename> uses only alphanumeric characters, hyphens, and "
            f"underscores (e.g., 'Sample01.fastq.gz', 'WW_SiteA-2024.fastq.gz')."
        )
    return valid_files


# ---------------------------------------------------------------------------
# Step 1 – Quality filtering and amplicon retrieval
# ---------------------------------------------------------------------------

def run_cutadapt(fastq_file, primer_csv, output_dir, num_cpus):
    """Trim ARTIC primer sequences and low-quality bases with Cutadapt.

    Reads primer sequences from *primer_csv* and trims both 5' adapters and
    bases with quality < 15 from both ends.  The output is written as a
    gzipped FASTQ file named ``<sample>_trimmed.fastq.gz``.

    Parameters
    ----------
    fastq_file : Path
        Input raw FASTQ file.
    primer_csv : Path
        CSV file with a ``cutadapt`` column listing primer sequences.
    output_dir : Path
        Directory for the trimmed output file.
    num_cpus : int
        Number of CPU cores to pass to Cutadapt (``-j``).

    Returns
    -------
    Path
        Path to the trimmed FASTQ file.
    """
    sample_name = fastq_file.stem
    if sample_name.endswith(".fastq"):
        sample_name = sample_name[:-6]
    trimmed = output_dir / f"{sample_name}_trimmed.fastq.gz"
    if trimmed.exists():
        return trimmed
    primers = ["-g"] + "|-g|".join(
        pd.read_csv(str(primer_csv)).cutadapt
    ).split("|")
    cmd = [
        "cutadapt", "-q", "15,15", "-j", str(num_cpus),
        "-o", str(trimmed), str(fastq_file),
    ] + primers
    subprocess.run(cmd, check=True, capture_output=True)
    return trimmed


def align_to_reference(trimmed_fastq, output_dir):
    """Align trimmed reads to the Wuhan-Hu-1 reference genome.

    Uses minimap2 with the Oxford Nanopore preset (``map-ont``).  The
    resulting SAM is converted to a coordinate-sorted, indexed BAM file.
    Intermediate SAM and unsorted BAM files are removed after sorting.

    Parameters
    ----------
    trimmed_fastq : Path
        Quality-trimmed FASTQ file produced by :func:`run_cutadapt`.
    output_dir : Path
        Directory for the sorted BAM output.

    Returns
    -------
    Path
        Path to the sorted and indexed BAM file.
    """
    sample_name = trimmed_fastq.stem.replace("_trimmed", "")
    if sample_name.endswith(".fastq"):
        sample_name = sample_name[:-6]
    sorted_bam = output_dir / f"{sample_name}.sorted.bam"
    if sorted_bam.exists():
        return sorted_bam
    sam_file = output_dir / f"{sample_name}.sam"
    bam_file = output_dir / f"{sample_name}.bam"
    with open(sam_file, "w") as fh:
        subprocess.run(
            ["minimap2", "-ax", "map-ont", str(REF_GENOME), str(trimmed_fastq)],
            stdout=fh, check=True,
        )
    with open(bam_file, "wb") as fh:
        subprocess.run(["samtools", "view", "-Sb", str(sam_file)], stdout=fh, check=True)
    subprocess.run(["samtools", "sort", "-o", str(sorted_bam), str(bam_file)], check=True)
    subprocess.run(["samtools", "index", str(sorted_bam)], check=True)
    sam_file.unlink()
    bam_file.unlink()
    return sorted_bam


def extract_amplicon_reads(sorted_bam, amplicon, region, output_dir):
    """Extract reads mapping to a specific amplicon region from a BAM file.

    SAMtools filters unmapped reads (``-F 4``) in the specified genomic
    *region* on NC_045512.2, then converts the filtered BAM back to FASTQ.
    Temporary BAM files are cleaned up after extraction.

    Parameters
    ----------
    sorted_bam : Path
        Coordinate-sorted BAM file from :func:`align_to_reference`.
    amplicon : str
        Amplicon identifier (e.g. ``"A72"``).
    region : str
        Genomic region string (e.g. ``"21774-21923"``).
    output_dir : Path
        Parent directory where an amplicon subdirectory will be created.

    Returns
    -------
    Path or None
        Path to the extracted FASTQ file, or ``None`` if no reads were found.
    """
    sample_name = sorted_bam.stem.replace(".sorted", "")
    if sample_name.endswith(".fastq"):
        sample_name = sample_name[:-6]
    amplicon_dir = output_dir / amplicon
    amplicon_dir.mkdir(exist_ok=True)
    fastq_out = amplicon_dir / f"{sample_name}.fastq"
    if fastq_out.exists():
        return fastq_out if fastq_out.stat().st_size > 0 else None
    filtered_bam = output_dir / f"{sample_name}.{amplicon}.filtered.bam"
    sorted_filtered = output_dir / f"{sample_name}.{amplicon}.sorted.filtered.bam"
    with open(filtered_bam, "wb") as fh:
        subprocess.run(
            ["samtools", "view", "-bh", "-F", "4", str(sorted_bam),
             f"NC_045512.2:{region}"],
            stdout=fh, check=True,
        )
    subprocess.run(
        ["samtools", "sort", "-o", str(sorted_filtered), str(filtered_bam)],
        check=True,
    )
    subprocess.run(["samtools", "index", str(sorted_filtered)], check=True)
    filtered_bam.unlink()
    if sorted_filtered.stat().st_size > 0:
        with open(fastq_out, "w") as fh:
            subprocess.run(["samtools", "fastq", str(sorted_filtered)], stdout=fh, check=True)
    sorted_filtered.unlink(missing_ok=True)
    Path(str(sorted_filtered) + ".bai").unlink(missing_ok=True)
    return fastq_out if fastq_out.exists() and fastq_out.stat().st_size > 0 else None


def step1_amplicon_retrieval(artic_name, config, num_cpus):
    """Step 1: Quality filtering and amplicon retrieval for one ARTIC version.

    For every valid FASTQ file in the ARTIC directory, this step:
      1. Trims primers and low-quality bases (Cutadapt).
      2. Aligns reads to the Wuhan-Hu-1 reference (minimap2).
      3. Extracts reads for each Spike amplicon region (SAMtools).

    Parameters
    ----------
    artic_name : str
        ARTIC version label (e.g. ``"ARTIC_3"``).
    config : dict
        Configuration dictionary from :data:`ARTIC_CONFIGS`.
    num_cpus : int
        Number of CPU cores for Cutadapt.

    Returns
    -------
    dict[str, list[Path]]
        Mapping of amplicon name to list of per-sample amplicon FASTQ files.
    """
    artic_dir = BASE_DIR / artic_name
    output_dir = artic_dir / "Aligned_reads"
    output_dir.mkdir(exist_ok=True)
    fastq_files = validate_fastq_filenames(artic_dir, artic_name)
    if not fastq_files:
        logger.warning(f"[{artic_name}] No valid FASTQ files found")
        return {}
    logger.info(f"[{artic_name}] Processing {len(fastq_files)} FASTQ files")
    sorted_bams = []
    for fq in fastq_files:
        trimmed = run_cutadapt(fq, config["primer_csv"], output_dir, num_cpus)
        bam = align_to_reference(trimmed, output_dir)
        sorted_bams.append(bam)
    amplicon_fastqs = {}
    for amplicon, region in config["positions"].items():
        amplicon_fastqs[amplicon] = []
        for bam in sorted_bams:
            fq = extract_amplicon_reads(bam, amplicon, region, output_dir)
            if fq is not None:
                amplicon_fastqs[amplicon].append(fq)
    return amplicon_fastqs

# ---------------------------------------------------------------------------
# Step 2 – Clustering and BLAST
# ---------------------------------------------------------------------------

def convert_fastq_to_fasta(fastq_file, output_dir):
    """Convert a FASTQ file to FASTA format using seqtk.

    Parameters
    ----------
    fastq_file : Path
        Input FASTQ file.
    output_dir : Path
        Directory for the output FASTA file.

    Returns
    -------
    Path or None
        Path to the FASTA file, or ``None`` if conversion failed.
    """
    sample = fastq_file.stem
    fasta_out = output_dir / f"{sample}.fasta"
    if fasta_out.exists():
        return fasta_out
    try:
        with open(fasta_out, "w") as fh:
            subprocess.run(
                ["seqtk", "seq", "-a", str(fastq_file)],
                stdout=fh, check=True,
            )
    except subprocess.CalledProcessError:
        logger.error(f"FASTQ-to-FASTA conversion failed for {fastq_file.name}")
        fasta_out.unlink(missing_ok=True)
        return None
    return fasta_out


def cluster_sequences(fasta_file, output_dir):
    """Cluster sequences at 100% identity using VSEARCH.

    Retains only sequences between 300 and 500 nucleotides.  Cluster
    consensus sequences are written with ``size=`` annotations for
    downstream abundance filtering.

    Parameters
    ----------
    fasta_file : Path
        Input FASTA file of amplicon reads.
    output_dir : Path
        Directory for the cluster consensus output.

    Returns
    -------
    Path
        Path to the cluster consensus FASTA file.
    """
    sample = fasta_file.stem
    cluster_out = output_dir / f"{sample}_clusters.fasta"
    if cluster_out.exists():
        return cluster_out
    subprocess.run([
        "vsearch", "--cluster_fast", str(fasta_file),
        "--id", "1", "--consout", str(cluster_out),
        "--minseqlength", "300", "--maxseqlength", "500",
        "--sizeout", "--sizein", "--clusterout_sort",
        "--strand", "both", "--threads", "16",
    ], check=True, capture_output=True)
    return cluster_out


def blast_against_reference(cluster_fasta, ref_fasta, output_dir):
    """Search cluster centroids against an amplicon nucleotide reference.

    Uses VSEARCH ``--usearch_global`` at a minimum 50% identity threshold.
    Pairwise alignment FASTA files and BLAST-6 tabular results are generated.

    Parameters
    ----------
    cluster_fasta : Path
        Cluster centroid FASTA from :func:`cluster_sequences`.
    ref_fasta : Path
        Amplicon nucleotide reference FASTA.
    output_dir : Path
        Directory for output files.

    Returns
    -------
    Path or None
        Path to the pairwise alignment FASTA, or ``None`` if no hits.
    """
    sample = cluster_fasta.stem
    results_txt = output_dir / f"{sample}_results.txt"
    pairwise = output_dir / f"{sample}_pairwise_alignments.fasta"
    matched = output_dir / f"{sample}_matched.fasta"
    notmatched = output_dir / f"{sample}_notmatched.fasta"
    alnout = output_dir / f"{sample}_alignments.txt"
    if results_txt.exists():
        return pairwise if pairwise.exists() and pairwise.stat().st_size > 0 else None
    subprocess.run([
        "vsearch", "--usearch_global", str(cluster_fasta),
        "--db", str(ref_fasta),
        "--id", "0.50", "--strand", "both",
        "--alnout", str(alnout), "--blast6out", str(results_txt),
        "--notmatched", str(notmatched), "--matched", str(matched),
        "--sizein", "--threads", "16",
        "--fastapairs", str(pairwise),
        "--top_hits_only", "--maxhits", "1",
    ], check=True, capture_output=True)
    if matched.exists() and matched.stat().st_size > 0 and results_txt.stat().st_size > 0:
        return pairwise
    return None


def step2_cluster_and_blast(amplicon_fastqs, artic_name, config):
    """Step 2: Sequence clustering and reference search for all amplicons.

    For each amplicon, converts per-sample FASTQ to FASTA, clusters at 100%
    identity, and searches cluster centroids against the amplicon nucleotide
    reference using VSEARCH.

    Parameters
    ----------
    amplicon_fastqs : dict[str, list[Path]]
        Output of :func:`step1_amplicon_retrieval`.
    artic_name : str
        ARTIC version label.
    config : dict
        Configuration dictionary from :data:`ARTIC_CONFIGS`.

    Returns
    -------
    dict[str, list[Path]]
        Mapping of amplicon name to list of pairwise alignment FASTA files.
    """
    artic_dir = BASE_DIR / artic_name
    ref_subdir = AMPLICON_REF_DIR / config["amplicon_ref_subdir"]
    pairwise_files = {}
    for amplicon, fastq_list in amplicon_fastqs.items():
        ref_fasta = ref_subdir / f"{amplicon}.fasta"
        if not ref_fasta.exists():
            logger.info(f"[{artic_name}/{amplicon}] No reference found, skipping")
            continue
        amplicon_dir = artic_dir / "Aligned_reads" / amplicon
        vsearch_dir = amplicon_dir / "vsearch_results"
        vsearch_dir.mkdir(exist_ok=True)
        pairwise_files[amplicon] = []
        for fq in fastq_list:
            fasta = convert_fastq_to_fasta(fq, vsearch_dir)
            clustered = cluster_sequences(fasta, vsearch_dir)
            if clustered.exists() and clustered.stat().st_size > 0:
                pw = blast_against_reference(clustered, ref_fasta, vsearch_dir)
                if pw and pw.exists():
                    pairwise_files[amplicon].append(pw)
    return pairwise_files

# ---------------------------------------------------------------------------
# Step 3 – Curate pairwise alignments
# ---------------------------------------------------------------------------

def get_reference_id(ref_fasta):
    """Return the sequence identifier of the first record in a FASTA file.

    Parameters
    ----------
    ref_fasta : Path
        Reference FASTA file.

    Returns
    -------
    str
        The FASTA header identifier (first whitespace-delimited token).
    """
    record = next(SeqIO.parse(str(ref_fasta), "fasta"))
    return record.id


def curate_pairwise_alignment(pairwise_fasta, ref_id, output_dir, sample_name):
    """Curate a VSEARCH pairwise alignment by removing the reference sequence.

    Reads the interleaved pairwise alignment FASTA, strips all records whose
    header matches *ref_id*, and removes gap characters (``-``) from the
    remaining query sequences to produce clean, ungapped consensus reads.

    Parameters
    ----------
    pairwise_fasta : Path
        Pairwise alignment FASTA from VSEARCH.
    ref_id : str
        Reference sequence identifier to remove.
    output_dir : Path
        Directory for the curated output.
    sample_name : str
        Sample identifier used to name the output file.

    Returns
    -------
    Path or None
        Path to the curated FASTA, or ``None`` if empty.
    """
    curated = output_dir / f"{sample_name}_curated.fasta"
    with open(pairwise_fasta, "r") as fh:
        lines = fh.readlines()
    with open(curated, "w") as fh:
        write_seq = True
        for line in lines:
            if line.startswith(">"):
                header_id = line.strip().lstrip(">").split()[0]
                if header_id == ref_id:
                    write_seq = False
                else:
                    write_seq = True
                    fh.write(line)
            else:
                if write_seq:
                    fh.write(line.replace("-", ""))
    if curated.stat().st_size > 0:
        return curated
    curated.unlink(missing_ok=True)
    return None


def step3_curate(pairwise_files, artic_name, config):
    """Step 3: Curate pairwise alignments for all amplicons.

    Removes reference sequences from VSEARCH pairwise alignment outputs and
    strips alignment gaps from query sequences.

    Parameters
    ----------
    pairwise_files : dict[str, list[Path]]
        Output of :func:`step2_cluster_and_blast`.
    artic_name : str
        ARTIC version label.
    config : dict
        Configuration dictionary from :data:`ARTIC_CONFIGS`.

    Returns
    -------
    dict[str, list[Path]]
        Mapping of amplicon name to list of curated FASTA files.
    """
    ref_subdir = AMPLICON_REF_DIR / config["amplicon_ref_subdir"]
    curated_files = {}
    for amplicon, pw_list in pairwise_files.items():
        ref_fasta = ref_subdir / f"{amplicon}.fasta"
        ref_id = get_reference_id(ref_fasta)
        curated_files[amplicon] = []
        for pw in pw_list:
            sample_name = pw.stem.replace("_clusters_pairwise_alignments", "")
            curated = curate_pairwise_alignment(pw, ref_id, pw.parent, sample_name)
            if curated:
                curated_files[amplicon].append(curated)
    return curated_files

# ---------------------------------------------------------------------------
# Step 4 – Filter by sequence count and translate to amino acids
# ---------------------------------------------------------------------------

def filter_by_seq_count(fasta_file, output_dir, sample_name, min_seqs=2):
    """Discard cluster centroids supported by fewer than *min_seqs* reads.

    Reads the ``seqs=`` annotation in each FASTA header and retains only
    records meeting the minimum abundance threshold.

    Parameters
    ----------
    fasta_file : Path
        Curated FASTA from Step 3.
    output_dir : Path
        Directory for the filtered output.
    sample_name : str
        Sample identifier for the output filename.
    min_seqs : int, default 2
        Minimum cluster size to retain.

    Returns
    -------
    Path or None
        Path to the filtered FASTA, or ``None`` if no records passed.
    """
    filtered = output_dir / f"{sample_name}_filtered.fasta"
    records = []
    for record in SeqIO.parse(str(fasta_file), "fasta"):
        match = re.search(r"seqs=(\d+)", record.id)
        if match and int(match.group(1)) >= min_seqs:
            records.append(record)
    if records:
        SeqIO.write(records, str(filtered), "fasta")
        return filtered
    return None


def translate_to_aa(fasta_file, output_dir, sample_name, is_last_amplicon=False):
    """Translate nucleotide sequences to amino acids in the best reading frame.

    Each sequence is translated in all three forward reading frames.  The
    longest ORF without internal stop codons is selected.  For the last
    amplicon of each ARTIC scheme a terminal stop codon is permitted (it
    represents the natural end of the Spike protein).  Sequences with
    excessive single-character (≥7) or dinucleotide (≥7) repeats are
    discarded as likely artefacts.

    Parameters
    ----------
    fasta_file : Path
        Filtered nucleotide FASTA from :func:`filter_by_seq_count`.
    output_dir : Path
        Directory for the translated output.
    sample_name : str
        Sample identifier used to name the output file.
    is_last_amplicon : bool, default False
        If ``True``, allow a terminal ``*`` in the amino acid sequence.

    Returns
    -------
    Path or None
        Path to the translated amino acid FASTA, or ``None`` if empty.
    """
    translated = output_dir / f"{sample_name}_translated.fasta"

    def valid_translation(seq, allow_terminal_stop):
        valid = []
        for frame in range(3):
            trimmed = seq[frame:]
            if len(trimmed) % 3 != 0:
                trimmed += "N" * (3 - len(trimmed) % 3)
            aa = str(Seq(trimmed).translate(to_stop=False)).rstrip("X")
            if allow_terminal_stop and aa.endswith("*"):
                aa = aa[:-1]
            if "*" not in aa:
                valid.append(aa)
        return max(valid, key=len) if valid else None

    def has_excessive_repeats(seq):
        if re.search(r"(.)\1{6,}", str(seq)):
            return True
        if re.search(r"(..)\1{6,}", str(seq)):
            return True
        return False

    records = []
    for record in SeqIO.parse(str(fasta_file), "fasta"):
        aa_seq = valid_translation(str(record.seq), allow_terminal_stop=is_last_amplicon)
        if aa_seq and not has_excessive_repeats(aa_seq):
            record.seq = Seq(aa_seq)
            records.append(record)
    if records:
        SeqIO.write(records, str(translated), "fasta")
        return translated
    return None


def step4_filter_and_translate(curated_files, config):
    """Step 4: Abundance filtering and amino acid translation for all amplicons.

    Discards low-abundance cluster centroids and translates the remaining
    nucleotide sequences to amino acids.

    Parameters
    ----------
    curated_files : dict[str, list[Path]]
        Output of :func:`step3_curate`.
    config : dict
        Configuration dictionary from :data:`ARTIC_CONFIGS`.

    Returns
    -------
    dict[str, list[Path]]
        Mapping of amplicon name to list of translated amino acid FASTA files.
    """
    last_amplicon = config["last_amplicon"]
    translated_files = {}
    for amplicon, curated_list in curated_files.items():
        is_last = amplicon == last_amplicon
        translated_files[amplicon] = []
        for curated in curated_list:
            sample_name = curated.stem.replace("_curated", "")
            filtered = filter_by_seq_count(curated, curated.parent, sample_name)
            if filtered:
                translated = translate_to_aa(
                    filtered, filtered.parent, sample_name, is_last_amplicon=is_last
                )
                if translated:
                    translated_files[amplicon].append(translated)
    return translated_files

# ---------------------------------------------------------------------------
# Step 5 – MMseqs2 alignment against spike database
# ---------------------------------------------------------------------------

def mmseqs_align(fasta_file, output_dir):
    """Align amino acid sequences against the Spike protein database with MMseqs2.

    Runs ``mmseqs easy-search`` with sensitivity 7.5, minimum 90% sequence
    identity, and minimum 90% query coverage.  Only the best hit per query
    is retained.  Sequences passing these thresholds are written to a FASTA
    file with alignment gaps removed.

    Parameters
    ----------
    fasta_file : Path
        Translated amino acid FASTA from Step 4.
    output_dir : Path
        Directory for alignment outputs (TSV results and aligned FASTA).

    Returns
    -------
    Path or None
        Path to the aligned amino acid FASTA, or ``None`` if no sequences
        passed the thresholds.
    """
    sample = fasta_file.stem
    temp_dir = output_dir / f"temp_{sample}"
    query_fasta = temp_dir / f"{sample}_query.fasta"
    result_tsv = output_dir / f"{sample}_result.tsv"
    aligned_fasta = output_dir / f"{sample}_aligned.fasta"
    if aligned_fasta.exists():
        return aligned_fasta
    temp_dir.mkdir(exist_ok=True)
    SeqIO.write(SeqIO.parse(str(fasta_file), "fasta"), str(query_fasta), "fasta")
    try:
        subprocess.run([
            MMSEQS_EXECUTABLE, "easy-search",
            str(query_fasta), str(SPIKE_DB), str(result_tsv), str(temp_dir),
            "--format-output",
            "query,target,qaln,taln,qstart,qend,tstart,tend,qlen,tlen,"
            "mismatch,gapopen,fident,nident,pident,qcov,tcov,evalue,bits,cigar,qseq,tseq",
            "-s", "7.5", "--alignment-mode", "0",
            "--min-seq-id", "0.9", "--cov-mode", "2", "-c", "0.9",
            "--max-seqs", "1",
        ], check=True, capture_output=True)
    except subprocess.CalledProcessError as exc:
        logger.error(f"MMseqs2 failed for {fasta_file.name}: {exc}")
        shutil.rmtree(temp_dir, ignore_errors=True)
        return None
    aligned_sequences = {}
    if result_tsv.exists():
        with open(result_tsv) as fh:
            for line in fh:
                fields = line.strip().split("\t")
                query, qaln = fields[0], fields[2]
                aligned_sequences[query] = qaln.replace("-", "")
    shutil.rmtree(temp_dir, ignore_errors=True)
    if aligned_sequences:
        with open(aligned_fasta, "w") as fh:
            for qid, seq in aligned_sequences.items():
                fh.write(f">{qid}\n{seq}\n")
        return aligned_fasta
    return None


def step5_mmseqs_alignment(translated_files, artic_name):
    """Step 5: MMseqs2 amino acid alignment for all amplicons.

    Aligns translated cluster centroids against the curated Spike protein
    database.  Only sequences meeting the identity and coverage thresholds
    are retained for downstream entropy analysis.

    Parameters
    ----------
    translated_files : dict[str, list[Path]]
        Output of :func:`step4_filter_and_translate`.
    artic_name : str
        ARTIC version label.

    Returns
    -------
    dict[str, list[Path]]
        Mapping of amplicon name to list of aligned amino acid FASTA files.
    """
    artic_dir = BASE_DIR / artic_name
    aligned_files = {}
    for amplicon, trans_list in translated_files.items():
        aa_dir = artic_dir / "Aligned_reads" / amplicon / "aa_alignment"
        aa_dir.mkdir(exist_ok=True)
        aligned_files[amplicon] = []
        for trans in trans_list:
            aligned = mmseqs_align(trans, aa_dir)
            if aligned:
                aligned_files[amplicon].append(aligned)
    return aligned_files

# ---------------------------------------------------------------------------
# Step 6 – Per-position entropy analysis
# ---------------------------------------------------------------------------

def parse_weight_from_header(header: str) -> int:
    """Extract the cluster-size weight from a FASTA header annotation.

    Looks for ``seqs=<N>`` first (VSEARCH consensus output), then falls back
    to ``size=<N>`` (VSEARCH centroid output).  Returns 1 if neither tag is
    found.

    Parameters
    ----------
    header : str
        Full FASTA description line.

    Returns
    -------
    int
        The cluster size weight.
    """
    match = re.search(r"seqs=(\d+)", header)
    if match:
        return int(match.group(1))
    match = re.search(r"size=(\d+)", header)
    if match:
        return int(match.group(1))
    return 1


def _clean_mafft_id(seq_id: str) -> str:
    """Return the first whitespace-delimited token of a MAFFT sequence ID."""
    return seq_id.split()[0]


def run_pairwise_mafft(
    ref_record: SeqRecord, query_record: SeqRecord, workdir_root: Path
) -> tuple:
    """Align a reference and a single query sequence using MAFFT.

    Creates a temporary FASTA with both sequences (labelled ``REF`` and
    ``QUERY``), runs MAFFT with ``--localpair --maxiterate 1000`` and custom
    gap penalties (OP=12, EP=3), and returns the aligned strings.

    Parameters
    ----------
    ref_record : SeqRecord
        Reference amino acid sequence.
    query_record : SeqRecord
        Query amino acid sequence.
    workdir_root : Path
        Parent directory for temporary files.

    Returns
    -------
    tuple[str, str]
        ``(aligned_ref, aligned_query)`` — gapped alignment strings.

    Raises
    ------
    RuntimeError
        If the MAFFT output does not contain both REF and QUERY records.
    """
    workdir_root = Path(workdir_root)
    with tempfile.TemporaryDirectory(prefix="mafft_pair_", dir=str(workdir_root)) as tmp_dir:
        input_fa = Path(tmp_dir) / "pair.fasta"
        ref_for_alignment = SeqRecord(ref_record.seq, id="REF", description="REF")
        query_for_alignment = SeqRecord(
            query_record.seq, id="QUERY", description=query_record.description
        )
        SeqIO.write([ref_for_alignment, query_for_alignment], input_fa, "fasta")
        cmd = [
            MAFFT_EXECUTABLE, "--localpair", "--maxiterate", "1000",
            "--op", "12", "--ep", "3", "--thread", "1", str(input_fa),
        ]
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        aligned_records = list(SeqIO.parse(io.StringIO(proc.stdout), "fasta"))
    ref_seq = None
    query_seq = None
    for rec in aligned_records:
        cleaned_id = _clean_mafft_id(rec.id)
        if cleaned_id == "REF" and ref_seq is None:
            ref_seq = str(rec.seq)
        elif cleaned_id == "QUERY" and query_seq is None:
            query_seq = str(rec.seq)
    if ref_seq is None or query_seq is None:
        found_ids = [rec.id for rec in aligned_records]
        raise RuntimeError(
            f"MAFFT output missing REF/QUERY sequences. Found IDs: {found_ids}"
        )
    return ref_seq, query_seq


def translate_amplicon_reference(nuc_path: Path) -> SeqRecord:
    """Translate an amplicon nucleotide reference to amino acids.

    The nucleotide sequence is trimmed to a multiple of three and translated
    in reading frame 0.  Any stop codons (``*``) are removed from the
    resulting amino acid sequence.

    Parameters
    ----------
    nuc_path : Path
        Path to the amplicon nucleotide FASTA file.

    Returns
    -------
    SeqRecord
        Translated amino acid record with ``_AA`` appended to the ID.

    Raises
    ------
    ValueError
        If the nucleotide sequence is too short for translation.
    """
    nuc_record = next(SeqIO.parse(str(nuc_path), "fasta"))
    nuc_seq = nuc_record.seq
    trimmed_len = len(nuc_seq) - (len(nuc_seq) % 3)
    if trimmed_len <= 0:
        raise ValueError(
            f"Amplicon reference {nuc_path} length is not sufficient for translation."
        )
    if trimmed_len != len(nuc_seq):
        nuc_seq = nuc_seq[:trimmed_len]
    aa_seq = nuc_seq.translate(to_stop=False)
    aa_seq = Seq(str(aa_seq).replace("*", ""))
    return SeqRecord(
        aa_seq, id=f"{nuc_record.id}_AA",
        description=f"{nuc_record.description} translated",
    )


def map_amplicon_positions_to_wuhan(
    wuhan_record: SeqRecord, amplicon_record: SeqRecord, workdir_root: Path
) -> list:
    """Map each amplicon amino acid position to a Wuhan-Hu-1 Spike position.

    Aligns the amplicon amino acid sequence to the full Wuhan-Hu-1 Spike
    protein using MAFFT, then walks the alignment columns to produce a
    one-to-one mapping from 0-based amplicon positions to 1-based Wuhan
    Spike positions.

    Parameters
    ----------
    wuhan_record : SeqRecord
        Wuhan-Hu-1 Spike amino acid reference.
    amplicon_record : SeqRecord
        Translated amplicon amino acid reference.
    workdir_root : Path
        Working directory for MAFFT temporary files.

    Returns
    -------
    list[int]
        List of Wuhan Spike positions (1-based) corresponding to each
        amplicon amino acid position.

    Raises
    ------
    ValueError
        If any amplicon position cannot be mapped to the Wuhan reference.
    """
    aligned_wuhan, aligned_amplicon = run_pairwise_mafft(
        wuhan_record, amplicon_record, workdir_root
    )
    mapping = [None] * len(amplicon_record.seq)
    wuhan_pos = 0
    amplicon_pos = 0
    for w_char, a_char in zip(aligned_wuhan, aligned_amplicon):
        if w_char not in GAP_CHARS:
            wuhan_pos += 1
        if a_char not in GAP_CHARS:
            amplicon_pos += 1
            if amplicon_pos > len(mapping):
                break
            if w_char in GAP_CHARS:
                raise ValueError(
                    "Amplicon amino acid aligns to a gap in Wuhan reference; "
                    "cannot determine position."
                )
            mapping[amplicon_pos - 1] = wuhan_pos
    if any(pos is None for pos in mapping):
        missing = ", ".join(
            str(idx + 1) for idx, pos in enumerate(mapping) if pos is None
        )
        raise ValueError(
            f"Unable to map all amplicon positions to Wuhan reference (missing: {missing})."
        )
    return [int(pos) for pos in mapping]


def persist_amplicon_reference(
    record: SeqRecord, output_dir: Path, amplicon_name: str
) -> Path:
    """Write the translated amplicon AA reference to disk.

    The file is consumed by worker processes during parallel entropy
    computation, avoiding the need to serialise SeqRecord objects.

    Parameters
    ----------
    record : SeqRecord
        Translated amplicon amino acid sequence.
    output_dir : Path
        Directory in which to write the FASTA file.
    amplicon_name : str
        Amplicon identifier (e.g. ``"A72"``) used in the filename.

    Returns
    -------
    Path
        Path to the written FASTA file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    target_path = output_dir / f"{amplicon_name}_AA_reference.fasta"
    SeqIO.write([record], target_path, "fasta")
    return target_path


def _trim_mplf_in_place(record: SeqRecord) -> None:
    """Remove all occurrences of the 4-AA insertion MPLF from a query AA record.

    This operates on the sequence as a string and removes every non-overlapping
    instance of the motif "MPLF" anywhere in the read.  Used to mitigate a
    problematic insertion at the beginning of ARTIC v5 amplicon A70 that causes
    alignment artifacts.
    """
    if not record.seq:
        return
    seq_str = str(record.seq)
    if "MPLF" not in seq_str:
        return
    record.seq = Seq(seq_str.replace("MPLF", ""))


def process_sample(
    fasta_path: Path,
    reference_record: SeqRecord,
    workdir_root: Path,
    position_labels: list,
    trim_mplf: bool = False,
) -> pd.DataFrame:
    """Build a per-read amino acid matrix for one sample FASTA.

    Each query read is individually aligned to the amplicon amino acid
    reference via MAFFT.  The alignment is walked column by column and the
    query residue (or ``"del"`` for a deletion) is recorded at each
    Wuhan Spike reference position.  Insertions relative to the reference
    are omitted.  Heuristics trim pathological head and tail alignments
    where a long run of deletions neighbours only a few residues at
    either end of the read.

    Parameters
    ----------
    fasta_path : Path
        Sample amino acid FASTA (cluster centroids with size annotations).
    reference_record : SeqRecord
        Translated amplicon amino acid reference.
    workdir_root : Path
        Working directory for MAFFT temporary files.
    position_labels : list[int]
        Wuhan Spike positions corresponding to each reference AA position.
    trim_mplf : bool, default False
        If ``True``, remove all ``MPLF`` motifs from each query before
        alignment (used for ARTIC v5 amplicon A70).

    Returns
    -------
    pd.DataFrame
        Wide matrix with columns ``ReadID``, ``Seqs``, and one column per
        Wuhan Spike position.
    """
    fasta_path = Path(fasta_path)
    ref_len = len(reference_record.seq)
    if ref_len != len(position_labels):
        raise ValueError("Reference length and position label length do not match.")
    records = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        if not rec.seq:
            continue
        records.append(rec)
    if not records:
        return pd.DataFrame(columns=["ReadID", "Seqs"] + position_labels)
    rows = []
    for rec in records:
        read_id = rec.id
        weight = parse_weight_from_header(rec.description or read_id)
        if trim_mplf:
            _trim_mplf_in_place(rec)
        aligned_ref, aligned_query = run_pairwise_mafft(
            reference_record, rec, workdir_root
        )
        if len(aligned_ref) != len(aligned_query):
            raise RuntimeError("Aligned reference and query have different lengths")
        first_q = None
        last_q = None
        for idx, ch in enumerate(aligned_query):
            if ch not in GAP_CHARS:
                first_q = idx
                break
        if first_q is None:
            continue
        for idx in range(len(aligned_query) - 1, -1, -1):
            if aligned_query[idx] not in GAP_CHARS:
                last_q = idx
                break
        row = {"ReadID": read_id, "Seqs": weight}
        events = []
        ref_pos = 0
        for col, (ref_char, q_char) in enumerate(zip(aligned_ref, aligned_query)):
            if ref_char in GAP_CHARS:
                continue
            ref_pos += 1
            if ref_pos > ref_len:
                break
            if col < first_q or col > last_q:
                continue
            symbol = "del" if q_char in GAP_CHARS else q_char
            column_label = position_labels[ref_pos - 1]
            row[column_label] = symbol
            events.append((ref_pos, column_label, symbol))

        # ----- Head-artifact heuristic -----
        # Mirror of the tail heuristic: if the alignment begins with very
        # few amino acids (≤ TAIL_NON_GAP_MAX) immediately followed by a
        # long run of deletions (≥ TAIL_GAP_MIN), the leading residues are
        # likely misplaced by MAFFT.  Clear the entire head region and
        # reinsert the amino acids just before the first reliably-aligned
        # position.
        head_non_gap_len = 0
        idx_h = 0
        while idx_h < len(events) and events[idx_h][2] != "del":
            head_non_gap_len += 1
            idx_h += 1
        head_gap_len = 0
        while idx_h < len(events) and events[idx_h][2] == "del":
            head_gap_len += 1
            idx_h += 1
        if (
            head_non_gap_len > 0
            and head_non_gap_len <= TAIL_NON_GAP_MAX
            and head_gap_len >= TAIL_GAP_MIN
        ):
            trim_total_h = head_non_gap_len + head_gap_len
            head_events = events[:trim_total_h]
            aa_events_h = [evt for evt in head_events if evt[2] != "del"]
            for _, col_label_trim, _ in head_events:
                row[col_label_trim] = ""
            if trim_total_h < len(events):
                first_kept_pos = events[trim_total_h][0]
                insert_pos = first_kept_pos - len(aa_events_h)
                for _, _, aa_symbol in aa_events_h:
                    if 1 <= insert_pos <= ref_len:
                        row[position_labels[insert_pos - 1]] = aa_symbol
                    insert_pos += 1

        # ----- Tail-artifact heuristic -----
        tail_non_gap_len = 0
        idx = len(events) - 1
        while idx >= 0 and events[idx][2] != "del":
            tail_non_gap_len += 1
            idx -= 1
        tail_gap_len = 0
        while idx >= 0 and events[idx][2] == "del":
            tail_gap_len += 1
            idx -= 1
        if (
            tail_non_gap_len > 0
            and tail_non_gap_len <= TAIL_NON_GAP_MAX
            and tail_gap_len >= TAIL_GAP_MIN
        ):
            trim_total = tail_non_gap_len + tail_gap_len
            tail_events = events[-trim_total:]
            aa_events = [evt for evt in tail_events if evt[2] != "del"]
            for _, col_label_trim, _ in tail_events:
                row[col_label_trim] = ""
            last_kept_pos = (
                events[-trim_total - 1][0] if len(events) > trim_total else 0
            )
            next_pos = last_kept_pos + 1
            for _, _, aa_symbol in aa_events:
                if next_pos > ref_len:
                    break
                row[position_labels[next_pos - 1]] = aa_symbol
                next_pos += 1
        for col_label in position_labels:
            row.setdefault(col_label, "")
        rows.append(row)
    columns = ["ReadID", "Seqs"] + position_labels
    return pd.DataFrame(rows, columns=columns)


def counts_from_matrix(df: pd.DataFrame, position_labels: list) -> dict:
    """Aggregate weighted amino acid counts per reference position.

    Iterates through the per-read matrix produced by :func:`process_sample`
    and accumulates symbol counts at each Wuhan Spike position, weighting
    each row by its ``Seqs`` (cluster size) value.  Empty cells and ``NaN``
    values are skipped.

    Parameters
    ----------
    df : pd.DataFrame
        Per-read amino acid matrix.
    position_labels : list[int]
        Wuhan Spike reference positions.

    Returns
    -------
    dict[int, Counter]
        Mapping of position to Counter of amino acid symbols.
    """
    counts = defaultdict(Counter)
    for label in position_labels:
        counts[label]
    if df.empty:
        return counts
    for _, row in df.iterrows():
        weight = row.get("Seqs", 1)
        if pd.isna(weight):
            weight = 1
        try:
            weight = int(weight)
        except (ValueError, TypeError):
            weight = 1
        for label in position_labels:
            if label not in row:
                continue
            val = row.get(label)
            if pd.isna(val) or val == "":
                continue
            counts[label][str(val)] += weight
    return counts


def counts_to_entropy_df(counts: dict, alpha: float = 0.01) -> pd.DataFrame:
    """Compute per-position entropy with minimal Dirichlet smoothing.

    Characteristics:
        * Only observed amino acids are smoothed.
        * No expansion to a fixed 20/21 AA alphabet.
        * Normalization by log2(min(21, total)) to reflect the maximum
          possible diversity given the read depth at this site.
    """
    columns = [
        "Position", "Entropy", "Most_Common", "Most_Common_Freq",
        "Coverage",
    ]
    positions = sorted(counts.keys())
    if not positions:
        return pd.DataFrame(columns=columns)
    rows = []
    for pos in positions:
        aa_counter = counts.get(pos, Counter())
        total = sum(aa_counter.values())
        if total == 0:
            rows.append({
                "Position": pos, "Entropy": "NA", "Most_Common": "NA",
                "Most_Common_Freq": "NA", "Coverage": "No",
            })
            continue
        observed_aas = [aa for aa in aa_counter.keys() if aa_counter[aa] > 0]
        K_obs = len(observed_aas)
        if K_obs == 0:
            rows.append({
                "Position": pos, "Entropy": "NA", "Most_Common": "NA",
                "Most_Common_Freq": "NA", "Coverage": "No",
            })
            continue
        denom = total + alpha * K_obs
        raw_entropy = 0.0
        for aa in observed_aas:
            p = (aa_counter[aa] + alpha) / denom
            raw_entropy -= p * math.log(p, 2)
        effective_K = min(21, total)
        if effective_K > 1:
            max_entropy = math.log(effective_K, 2)
            entropy = raw_entropy / max_entropy
        else:
            entropy = 0.0
        most_common_aa, most_common_count = aa_counter.most_common(1)[0]
        rows.append({
            "Position": pos, "Entropy": entropy,
            "Most_Common": most_common_aa,
            "Most_Common_Freq": most_common_count / total,
            "Coverage": "Yes",
        })
    return pd.DataFrame(rows)


def write_sample_excel(df: pd.DataFrame, out_path: Path) -> None:
    """Write a per-read amino acid matrix DataFrame to an Excel workbook.

    Parameters
    ----------
    df : pd.DataFrame
        Per-read matrix from :func:`process_sample`.
    out_path : Path
        Destination Excel file path.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="per_read_mapping")


def combine_csv_to_excel(output_dir: Path):
    """Aggregate per-sample entropy CSVs into one combined Excel workbook.

    Reads all ``*_entropy.csv`` files in *output_dir*, extracts the sample
    identifier from the filename, and produces a wide table with one row per
    sample and one column (``Pos_<N>``) per Wuhan Spike position.  The
    output is saved as ``combined_entropy_by_sample.xlsx``.

    Parameters
    ----------
    output_dir : Path
        Directory containing per-sample entropy CSV files.
    """
    csv_files = list(output_dir.glob("*_entropy.csv"))
    if not csv_files:
        return
    combined_data = {}
    all_positions = set()
    for csv_file in csv_files:
        identifier_match = re.match(r"(.+?)_translated", csv_file.stem)
        if not identifier_match:
            continue
        identifier = identifier_match.group(1)
        df = pd.read_csv(csv_file)
        entropy_dict = dict(zip(df["Position"], df["Entropy"]))
        combined_data[identifier] = entropy_dict
        all_positions.update(df["Position"])
    if not combined_data:
        return
    all_positions = sorted(all_positions)
    rows = []
    for identifier, entropy_dict in combined_data.items():
        row = {"Sample_ID": identifier}
        for pos in all_positions:
            row[f"Pos_{pos}"] = entropy_dict.get(pos, "NA")
        rows.append(row)
    final_df = pd.DataFrame(rows)
    excel_file = output_dir / "combined_entropy_by_sample.xlsx"
    final_df.to_excel(excel_file, index=False)
    logger.info(f"Combined entropy saved to {excel_file}")


def _process_sample_entry(
    sample_fasta: Path,
    amplicon_reference_path: Path,
    position_labels: list,
    output_dir: Path,
    trim_mplf: bool = False,
) -> str:
    """Worker entry point: process one sample and write outputs.

    Loads the amplicon reference from disk, builds the per-read amino acid
    matrix, computes Dirichlet-smoothed entropy, and writes the matrix
    Excel file and the entropy CSV.  Designed to run inside a
    :class:`~concurrent.futures.ProcessPoolExecutor`.

    Parameters
    ----------
    sample_fasta : Path
        Amino acid FASTA for one sample.
    amplicon_reference_path : Path
        On-disk path to the translated amplicon AA reference FASTA.
    position_labels : list[int]
        Wuhan Spike positions for this amplicon.
    output_dir : Path
        Directory for output Excel and CSV files.
    trim_mplf : bool, default False
        If ``True``, trim ``MPLF`` motifs before alignment.

    Returns
    -------
    str
        Human-readable status message.
    """
    reference_record = next(SeqIO.parse(str(amplicon_reference_path), "fasta"))
    workdir_root = amplicon_reference_path.parent
    df = process_sample(
        sample_fasta, reference_record, workdir_root, position_labels,
        trim_mplf=trim_mplf,
    )
    out_xlsx = output_dir / f"{Path(sample_fasta).stem}_per_read_matrix.xlsx"
    write_sample_excel(df, out_xlsx)
    counts = counts_from_matrix(df, position_labels)
    entropy_df = counts_to_entropy_df(counts)
    entropy_csv = output_dir / f"{Path(sample_fasta).stem}_entropy.csv"
    entropy_df.to_csv(entropy_csv, index=False)
    return f"{sample_fasta.name} -> {out_xlsx.name} & {entropy_csv.name}"


def step6_entropy(aligned_files, artic_name, config):
    """Step 6: Per-position amino acid entropy analysis for all amplicons.

    For each amplicon the function:
      1. Translates the amplicon nucleotide reference to amino acids.
      2. Maps amplicon positions to Wuhan-Hu-1 Spike coordinates.
      3. Aligns each query read to the amplicon reference (MAFFT).
      4. Builds per-read AA matrices and computes Dirichlet-smoothed entropy.
      5. Writes per-sample Excel/CSV files and a combined entropy summary.

    MPLF insertion trimming is automatically enabled for ARTIC v5 amplicon
    A70 to remove a known artefactual 4-AA insertion.

    Parameters
    ----------
    aligned_files : dict[str, list[Path]]
        Output of :func:`step5_mmseqs_alignment`.
    artic_name : str
        ARTIC version label.
    config : dict
        Configuration dictionary from :data:`ARTIC_CONFIGS`.
    """
    artic_dir = BASE_DIR / artic_name
    ref_subdir = AMPLICON_REF_DIR / config["amplicon_ref_subdir"]
    wuhan_record = next(SeqIO.parse(str(WUHAN_AA_SPIKE), "fasta"))
    for amplicon, sample_fastas in aligned_files.items():
        if not sample_fastas:
            continue
        output_dir = artic_dir / "Aligned_reads" / amplicon / "entropy"
        output_dir.mkdir(exist_ok=True)
        amplicon_ref_nuc = ref_subdir / f"{amplicon}.fasta"
        amplicon_aa_record = translate_amplicon_reference(amplicon_ref_nuc)
        amplicon_ref_path = persist_amplicon_reference(
            amplicon_aa_record, output_dir, amplicon
        )
        position_labels = map_amplicon_positions_to_wuhan(
            wuhan_record, amplicon_aa_record, workdir_root=output_dir
        )
        copied = []
        for fa in sample_fastas:
            target = output_dir / fa.name
            if not target.exists():
                shutil.copy2(fa, target)
            copied.append(target)
        wuhan_copy = output_dir / "wuhan_AA_Spike.fasta"
        if not wuhan_copy.exists():
            shutil.copy2(WUHAN_AA_SPIKE, wuhan_copy)
        # Enable MPLF trimming only for ARTIC_5 amplicon A70 to remove
        # a 4-AA insertion at the amplicon start that causes alignment artifacts.
        trim_mplf = (artic_name == "ARTIC_5" and amplicon == "A70")
        if trim_mplf:
            logger.info(
                f"[{artic_name}/{amplicon}] MPLF trimming enabled for this amplicon"
            )
        logger.info(
            f"[{artic_name}/{amplicon}] Running entropy on {len(copied)} samples"
        )
        max_workers = max(1, os.cpu_count() or 1)
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    _process_sample_entry,
                    fa, amplicon_ref_path, position_labels, output_dir,
                    trim_mplf,
                ): fa
                for fa in copied
            }
            for future in as_completed(futures):
                fa = futures[future]
                try:
                    msg = future.result()
                    logger.info(f"[{artic_name}/{amplicon}] {msg}")
                except Exception as exc:
                    logger.error(f"[{artic_name}/{amplicon}] {fa.name}: {exc}")
        combine_csv_to_excel(output_dir)

# ---------------------------------------------------------------------------
# Step 7 – Cross-amplicon entropy summary
# ---------------------------------------------------------------------------

def _rename_sample(sample: str) -> str:
    """Extract a short sample label from the original sample identifier.

    Recognised patterns (tried in order):
      1. ``YYMMDD-EDAR_NN-X``  → returns the numeric portion ``NN``
      2. ``YYMMDD_EDAR_NN``    → returns ``NN``
      3. ``YYMMDD-EDAR_NN``    → returns ``NN``
      4. ``<word>-YYYY-MM-DD``  → returns the leading word

    If none of these patterns match, the original identifier is returned
    unchanged.
    """
    if re.match(r"\d{6}-EDAR_\d{2,}-\w", sample):
        return sample.split("_")[1]
    if re.match(r"\d{6}_EDAR_\d{2,}", sample):
        return sample.split("_")[-1]
    if re.match(r"\d{6}-EDAR_\d{2,}", sample):
        return sample.split("_")[-1]
    if re.match(r"\w+-\d{4}-\d{2}-\d{2}", sample):
        return sample.split("-")[0]
    return sample


def _parse_date(filename: str) -> str:
    """Try to extract a date string (YYYY-MM-DD) from a sample identifier.

    Recognised patterns:
      1. ``D<word>-YYYY-MM-DD``  → returns ``YYYY-MM-DD``
      2. Any 6-digit group ``YYMMDD`` → converted to ``YYYY-MM-DD``

    Returns an empty string if no date can be determined.
    """
    if re.match(r"D\w+-\d{4}-\d{2}-\d{2}", filename):
        date_match = re.search(r"D\w+-(\d{4}-\d{2}-\d{2})", filename)
        if date_match:
            return date_match.group(1)
    date_match = re.search(r"(\d{6})", filename)
    if date_match:
        try:
            return datetime.strptime(date_match.group(1), "%y%m%d").strftime("%Y-%m-%d")
        except ValueError:
            pass
    return ""


def step7_generate_summary():
    """Merge per-amplicon combined entropy Excel files into one summary workbook.

    Produces ``per_position_entropy.xlsx`` in the pipeline root directory with
    one sheet per ARTIC version.  Each row is a sample; columns are Wuhan Spike
    amino acid positions.

    For ARTIC v5, amplicons A69 and A70 share overlapping Wuhan Spike
    positions.  The first 16 positions of A69 are given priority: entropy
    values from A70 (or any later amplicon) at those positions are discarded
    so that A69's values are preserved.  Non-overlapping positions from A70
    are retained normally.
    """
    artic_data: dict[str, dict] = {"3": {}, "4": {}, "5": {}}

    for artic_name in ["ARTIC_3", "ARTIC_4", "ARTIC_5"]:
        version = artic_name.split("_")[1]
        aligned_reads_dir = BASE_DIR / artic_name / "Aligned_reads"
        if not aligned_reads_dir.is_dir():
            continue

        # For ARTIC_5: the first 16 positions of A69 take priority over
        # overlapping positions in A70.  Populated when A69 is processed.
        a69_first16: set[int] = set()

        for amplicon_dir in sorted(aligned_reads_dir.iterdir()):
            if not amplicon_dir.is_dir() or not amplicon_dir.name.startswith("A"):
                continue
            amplicon_name = amplicon_dir.name
            entropy_dir = amplicon_dir / "entropy"
            if not entropy_dir.is_dir():
                continue
            combined_xlsx = entropy_dir / "combined_entropy_by_sample.xlsx"
            if not combined_xlsx.exists():
                continue

            try:
                df = pd.read_excel(combined_xlsx, engine="openpyxl")
            except Exception as exc:
                logger.warning(
                    f"[Summary] Could not read {combined_xlsx}: {exc}"
                )
                continue

            df = df.replace("NA", "")

            # Identify A69's first 16 Wuhan Spike positions for ARTIC_5
            if version == "5" and amplicon_name == "A69":
                all_pos_nums = sorted(
                    int(str(c).split("_")[1])
                    for c in df.columns
                    if str(c).startswith("Pos_")
                )
                a69_first16 = set(all_pos_nums[:16])

            for _, row in df.iterrows():
                try:
                    original_id = str(row["Sample_ID"])
                except KeyError:
                    continue

                if original_id not in artic_data[version]:
                    artic_data[version][original_id] = {
                        "Sample_ID": original_id,
                        "Sample": _rename_sample(original_id),
                        "Date": _parse_date(original_id),
                    }

                pos_columns = [c for c in df.columns if str(c).startswith("Pos_")]
                for col in pos_columns:
                    pos_num = int(str(col).split("_")[1])
                    # For ARTIC_5: skip positions already claimed by A69
                    if (
                        version == "5"
                        and amplicon_name != "A69"
                        and pos_num in a69_first16
                    ):
                        continue
                    value = row[col]
                    if pd.notna(value) and value != "":
                        artic_data[version][original_id][
                            f"Position_{pos_num}"
                        ] = value

    output_path = BASE_DIR / "per_position_entropy.xlsx"
    any_data = False
    try:
        with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
            for version in ["3", "4", "5"]:
                samples = artic_data[version]
                if not samples:
                    continue
                any_data = True
                df = pd.DataFrame.from_dict(samples, orient="index")
                pos_cols = sorted(
                    [c for c in df.columns if c.startswith("Position_")],
                    key=lambda x: int(x.split("_")[1]),
                )
                meta_cols = ["Sample_ID", "Sample", "Date"]
                df = df[meta_cols + pos_cols]
                df[pos_cols] = df[pos_cols].where(df[pos_cols].notna(), other="")
                df.to_excel(writer, sheet_name=f"ARTIC_{version}", index=False)
    except Exception as exc:
        logger.error(f"[Summary] Failed to write {output_path}: {exc}")
        return

    if any_data:
        logger.info(f"Summary entropy workbook saved to {output_path}")
    else:
        logger.warning("[Summary] No entropy data found; summary file not created.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Run the full integrated SARS-CoV-2 Spike amplicon pipeline.

    Iterates over all configured ARTIC primer schemes (v3, v4.1, v5.3.2)
    and for each one sequentially executes Steps 1 through 6 (quality
    filtering, clustering, curation, translation, MMseqs2 alignment, and
    entropy analysis).  After all versions are processed, Step 7 merges
    per-amplicon entropy results into a single summary workbook.

    Command-line arguments
    ----------------------
    -c / --cpus : int, optional
        Number of CPU cores (default: all available).
    """
    parser = argparse.ArgumentParser(
        description="Integrated SARS-CoV-2 Spike Amplicon Pipeline"
    )
    parser.add_argument(
        "-c", "--cpus", type=int, default=cpu_count(),
        help="Number of CPU cores to use (default: all)",
    )
    args = parser.parse_args()
    num_cpus = args.cpus

    for artic_name, config in ARTIC_CONFIGS.items():
        logger.info(f"{'='*60}")
        logger.info(f"Processing {artic_name}")
        logger.info(f"{'='*60}")

        artic_dir = BASE_DIR / artic_name
        if not artic_dir.is_dir():
            logger.info(f"[{artic_name}] Directory not found, skipping")
            continue
        if not list(artic_dir.glob("*.fastq.gz")):
            logger.info(f"[{artic_name}] No .fastq.gz files found, skipping")
            continue

        amplicon_fastqs = step1_amplicon_retrieval(artic_name, config, num_cpus)

        pairwise_files = step2_cluster_and_blast(amplicon_fastqs, artic_name, config)

        curated_files = step3_curate(pairwise_files, artic_name, config)

        translated_files = step4_filter_and_translate(curated_files, config)

        aligned_files = step5_mmseqs_alignment(translated_files, artic_name)

        step6_entropy(aligned_files, artic_name, config)

        logger.info(f"[{artic_name}] Finished")

    logger.info("Generating cross-amplicon entropy summary")
    step7_generate_summary()

    logger.info("Pipeline complete")


if __name__ == "__main__":
    main()
