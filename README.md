# Integrated Spike Entropy Pipeline

An end-to-end pipeline for processing SARS-CoV-2 Oxford Nanopore amplicon sequencing data targeting the Spike protein. The pipeline processes raw FASTQ reads from ARTIC primer schemes (v3, v4.1, and v5.3.2), performs quality filtering, amplicon extraction, sequence clustering, amino acid translation, database alignment, and per-position Shannon entropy calculation. Final outputs are per-sample Excel files containing amino acid matrices mapped to Wuhan-Hu-1 Spike reference positions and per-position entropy values.

---

## Table of Contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Input Data Structure](#input-data-structure)
4. [Reference Files](#reference-files)
5. [Usage](#usage)
6. [Pipeline Steps](#pipeline-steps)
7. [Output Structure](#output-structure)
8. [Parameters](#parameters)

---

## Requirements

The pipeline requires a Linux environment with [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) installed. All software dependencies are specified in the provided `environment.yml` file. The pipeline relies on the following tools and libraries:

### Command-line bioinformatics tools

| Tool     | Version  | Purpose                                                              |
| -------- | -------- | -------------------------------------------------------------------- |
| Cutadapt | 4.8      | Primer trimming and quality filtering                                |
| minimap2 | 2.28     | Long-read alignment to the reference genome                          |
| SAMtools | 1.21     | BAM file manipulation and amplicon region extraction                 |
| VSEARCH  | 2.28.1   | FASTQ-to-FASTA conversion, sequence clustering, and reference search |
| MMseqs2  | 15.6f452 | Amino acid sequence alignment against the Spike database             |
| MAFFT    | 7.526    | Pairwise amino acid alignment for positional entropy analysis        |

### Python libraries

| Library    | Version | Purpose                                                     |
| ---------- | ------- | ----------------------------------------------------------- |
| Python     | 3.10.14 | Runtime                                                     |
| pandas     | 2.2.2   | Data manipulation and CSV/Excel I/O                         |
| Biopython  | 1.85    | FASTA parsing, sequence translation, and SeqRecord handling |
| XlsxWriter | 3.2.5   | Excel file generation                                       |
| openpyxl   | 3.1.2   | Excel file support                                          |

---

## Installation

### 1. Install Conda

If Conda is not already installed, download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Miniforge](https://github.com/conda-forge/miniforge) (which includes Mamba by default):

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

Follow the prompts and restart the terminal session after installation.

### 2. Create the environment

Navigate to the directory containing the pipeline files and create the Conda environment from the provided `environment.yml`:

```bash
cd /path/to/pipeline
conda env create -f environment.yml
```

Alternatively, if Mamba is available (faster dependency resolution):

```bash
mamba env create -f environment.yml
```

This creates a Conda environment named `sars_spike_pipeline` with all required tools and Python packages at their specified versions.

### 3. Activate the environment

```bash
conda activate sars_spike_pipeline
```

### 4. Verify the installation

Confirm that all required tools are accessible:

```bash
cutadapt --version
minimap2 --version
samtools --version
vsearch --version
mmseqs version
mafft --version
python -c "import pandas; import Bio; import xlsxwriter; import openpyxl; print('All Python packages OK')"
```

---

## Input Data Structure

The pipeline expects the following directory layout in the same directory as the script:

```
pipeline_directory/
├── integrated_sars_pipeline.py
├── environment.yml
├── WG_wuhan.fa
├── wuhan_AA_Spike.fasta
├── sars_spike_database.fasta
├── Primer_schemes/
│   ├── primers_ARTICv.3.csv
│   ├── primers_ARTICv.4.1.csv
│   └── primers_ARTICv.5.3.2.csv
├── Amplicon_references/
│   ├── ARTIC3/
│   │   ├── A71.fasta
│   │   ├── A72.fasta
│   │   └── ...
│   ├── ARTIC4/
│   │   ├── A72.fasta
│   │   └── ...
│   └── ARTIC5/
│       ├── A70.fasta
│       └── ...
├── ARTIC_3/
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── ...
├── ARTIC_4/
│   ├── sample1.fastq.gz
│   └── ...
└── ARTIC_5/
    ├── sample1.fastq.gz
    └── ...
```

Place raw FASTQ files (gzipped, `.fastq.gz`) from each ARTIC primer scheme into its corresponding directory (`ARTIC_3/`, `ARTIC_4/`, or `ARTIC_5/`). Each file represents one sequencing sample.

> **Important — Filename requirements.** Input FASTQ filenames **must** follow the pattern `<samplename>.fastq.gz` where `<samplename>` contains only alphanumeric characters, hyphens (`-`), and underscores (`_`). Spaces, dots, or other special characters in the sample name portion are not supported. Additionally, sample names must not end with reserved pipeline suffixes (`_trimmed`, `_clusters`, `_curated`, `_filtered`, `_translated`, `_aligned`). Files that do not conform to this naming convention will be skipped with a warning message. Examples of valid filenames:
>
> - `Sample01.fastq.gz`
> - `Patient_A-2024.fastq.gz`
> - `barcode12.fastq.gz`
>
> Examples of **invalid** filenames:
>
> - `Sample 01.fastq.gz` (contains a space)
> - `Sample.01.fastq.gz` (contains a dot in the sample name)
> - `my_sample_trimmed.fastq.gz` (ends with reserved suffix `_trimmed`)

---

## Reference Files

| File                          | Description                                                                                       |
| ----------------------------- | ------------------------------------------------------------------------------------------------- |
| `WG_wuhan.fa`               | SARS-CoV-2 Wuhan-Hu-1 whole-genome reference (NC_045512.2) used for read alignment                |
| `wuhan_AA_Spike.fasta`      | Wuhan-Hu-1 Spike protein amino acid sequence used for positional mapping in the entropy step      |
| `sars_spike_database.fasta` | Curated SARS-CoV-2 Spike amino acid database for MMseqs2 variant alignment                        |
| `Primer_schemes/*.csv`      | Primer sequences for each ARTIC version in Cutadapt-compatible format                             |
| `Amplicon_references/`      | Per-amplicon nucleotide reference sequences for each ARTIC scheme, used as VSEARCH search targets |

---

## Usage

Activate the environment and run the pipeline from the directory containing the script:

```bash
conda activate sars_spike_pipeline
python integrated_sars_pipeline.py
```

To specify the number of CPU cores:

```bash
python integrated_sars_pipeline.py -c 8
```

By default, all available CPU cores are used.

---

## Pipeline Steps

The pipeline processes each ARTIC primer scheme (`ARTIC_3`, `ARTIC_4`, `ARTIC_5`) sequentially through six stages:

### Step 1 — Quality filtering and amplicon retrieval

For each sample FASTQ file:

1. **Primer trimming and quality filtering**: Cutadapt removes ARTIC primer sequences (loaded from the corresponding primer CSV) and trims bases with quality scores below 15 from both ends.
2. **Reference alignment**: Trimmed reads are aligned to the Wuhan-Hu-1 genome (`WG_wuhan.fa`) using minimap2 with Oxford Nanopore preset (`map-ont`).
3. **Amplicon extraction**: SAMtools extracts reads mapping to each amplicon region (ARTIC-version-specific genomic coordinates covering the Spike gene) and converts them back to FASTQ.

Each ARTIC version defines different amplicon boundaries. ARTIC v3 covers amplicons A71–A84; ARTIC v4 covers amplicons A72–A84; ARTIC v5 covers amplicons A70–A82.

### Step 2 — Sequence clustering and reference search

For each amplicon FASTQ:

1. **FASTQ-to-FASTA conversion**: VSEARCH converts reads from FASTQ to FASTA format.
2. **Clustering**: VSEARCH clusters sequences at 100% identity (`--id 1`), retaining only sequences between 300 and 500 nucleotides. Cluster consensus sequences are generated with size annotations.
3. **Reference search**: Cluster consensus sequences are searched against the corresponding amplicon nucleotide reference using VSEARCH global alignment (`--usearch_global`) at a minimum 50% identity threshold. Pairwise alignment files are generated for downstream curation.

### Step 3 — Alignment curation

Pairwise alignment FASTA files from Step 2 contain interleaved reference and query sequences. This step:

1. Removes reference sequences from each pairwise alignment (identified by matching the reference header from the amplicon FASTA).
2. Strips alignment gap characters (`-`) from the remaining query sequences, producing clean, ungapped consensus sequences.

### Step 4 — Sequence filtering and amino acid translation

1. **Abundance filtering**: Sequences with a cluster size (`seqs=`) below 2 are discarded.
2. **Translation**: Nucleotide sequences are translated to amino acids in all three forward reading frames. The longest translation without internal stop codons is retained. For the last amplicon of each ARTIC scheme (A84 for ARTIC v3/v4; A82 for ARTIC v5), a terminal stop codon is permitted as it represents the natural end of the Spike protein.
3. **Repeat filtering**: Translated sequences containing excessive single-character (≥7) or dinucleotide (≥7) repeats are discarded.

### Step 5 — MMseqs2 amino acid alignment

Translated amino acid sequences are aligned against the Spike protein database (`sars_spike_database.fasta`) using MMseqs2 `easy-search` with the following thresholds:

- Sensitivity: 7.5
- Minimum sequence identity: 90%
- Minimum query coverage: 90%
- Best hit only per query

Only sequences passing these thresholds are retained for entropy analysis.

### Step 6 — Per-position entropy analysis

For each amplicon and sample:

1. **Amplicon reference translation**: The amplicon nucleotide reference is translated to amino acids.
2. **Position mapping**: The amplicon amino acid sequence is aligned to the Wuhan-Hu-1 Spike amino acid reference (`wuhan_AA_Spike.fasta`) using MAFFT (`--localpair --maxiterate 1000`) to establish a mapping from amplicon positions to Wuhan Spike positions.
3. **Per-read alignment**: Each query amino acid sequence is individually aligned to the amplicon amino acid reference using MAFFT.
4. **Matrix construction**: A wide matrix is built with one row per read and one column per Wuhan Spike reference position. Cells contain the query amino acid at that position, `del` for deletions, or empty for positions outside the read coverage.
5. **Entropy calculation**: Weighted Shannon entropy (normalized by log₂ of the total weighted read count) is computed at each reference position.
6. **Output generation**: Per-sample Excel files with the full amino acid matrix and per-sample entropy CSV files are generated. A combined Excel file aggregates entropy values across all samples for each amplicon.

---

## Output Structure

After a complete run, the following output directories are created within each ARTIC folder:

```
ARTIC_X/
└── Aligned_reads/
    ├── sample.sorted.bam
    ├── sample_trimmed.fastq.gz
    ├── A7X/
    │   ├── sample.fastq
    │   ├── vsearch_results/
    │   │   ├── sample.fasta
    │   │   ├── sample_clusters.fasta
    │   │   ├── sample_clusters_pairwise_alignments.fasta
    │   │   ├── sample_clusters_matched.fasta
    │   │   ├── sample_curated.fasta
    │   │   ├── sample_filtered.fasta
    │   │   └── sample_translated.fasta
    │   ├── aa_alignment/
    │   │   ├── sample_translated_aligned.fasta
    │   │   └── sample_translated_result.tsv
    │   └── entropy/
    │       ├── A7X_AA_reference.fasta
    │       ├── wuhan_AA_Spike.fasta
    │       ├── sample_translated_aligned_per_read_matrix.xlsx
    │       ├── sample_translated_aligned_entropy.csv
    │       └── combined_entropy_by_sample.xlsx
    └── ...
```

### Key output files

| File                                | Description                                                                                                                  |
| ----------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `*_per_read_matrix.xlsx`          | Per-sample Excel file with one row per cluster centroid and columns for each Wuhan Spike amino acid position                 |
| `*_entropy.csv`                   | Per-sample CSV file with Shannon entropy, most common amino acid, and coverage at each position                              |
| `combined_entropy_by_sample.xlsx` | Aggregated entropy values across all samples for a given amplicon, with rows as samples and columns as Wuhan Spike positions |

---

## Parameters

| Argument           | Default             | Description                                              |
| ------------------ | ------------------- | -------------------------------------------------------- |
| `-c`, `--cpus` | All available cores | Number of CPU cores for Cutadapt and parallel processing |

Internal thresholds are defined in the script and can be modified if needed:

| Parameter                  | Value       | Context                                                 |
| -------------------------- | ----------- | ------------------------------------------------------- |
| Quality trimming threshold | 15          | Cutadapt base quality cutoff (both ends)                |
| Clustering identity        | 1.0 (100%)  | VSEARCH clustering threshold                            |
| Min/Max sequence length    | 300–500 nt | VSEARCH clustering length filter                        |
| VSEARCH search identity    | 0.50 (50%)  | Minimum identity for amplicon reference matching        |
| Minimum cluster size       | 2           | Sequences with fewer reads are discarded                |
| MMseqs2 sequence identity  | 0.9 (90%)   | Minimum identity for Spike database alignment           |
| MMseqs2 query coverage     | 0.9 (90%)   | Minimum coverage for Spike database alignment           |
| MAFFT gap penalties        | OP=12, EP=3 | Gap open and extension penalties for entropy alignments |
