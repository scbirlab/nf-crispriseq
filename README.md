# 🧮 CRISPRi-seq pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-crispriseq/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**scbirlab/nf-crispriseq** is a Nextflow pipeline that takes raw FASTQ files from a CRISPRi pooled screen and returns annotated per-guide read counts — and optionally per-guide fitness scores from a time course.

**Table of contents**

- [Quick start](#quick-start)
- [Worked example](#worked-example)
- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Fitness calculation](#fitness-calculation)
- [Running on a cluster](#running-on-a-cluster)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

---

## Quick start

> 🧬 **New here?** This pipeline takes raw sequencing files from a CRISPRi pooled screen and counts how many reads in each sample match each guide RNA. With time-course data, it can also score which gene knockdowns cause growth defects.

**1. Install Nextflow** (skip if already installed):

```bash
conda install -c bioconda nextflow
```

> On an HPC cluster (including the Crick), load modules instead: `module load Singularity Nextflow`

**2. Create a [sample sheet](#sample-sheet)** and, optionally, a [`nextflow.config`](#inputs) in your working directory.

**3. Run:**

```bash
nextflow run scbirlab/nf-crispriseq -r dev
```

First run installs dependencies from `environment.yml` (~10 min). To resume an interrupted run:

```bash
nextflow run scbirlab/nf-crispriseq -r dev -resume
```

---

## Worked example

This example processes a _Mycobacterium tuberculosis_ CRISPRi-seq dataset pulled directly from SRA.

### Directory layout

```
my-analysis/
├── nextflow.config
└── inputs/
    ├── sample-sheet.csv
    └── guides.fasta
```

### `inputs/guides.fasta`

One entry per guide RNA in the library:

```
>TargetA_guide1
TCGACTGAGCTGAAAGAAT
>TargetA_guide2
GTTTAAGAGCTATGCTGGT
>empty-vector
AACGTTAGCTAGCGATCGA
```

Alternatively, provide a CSV with `Name` and `guide_sequence` columns — see [Guide file formats](#guide-file-formats).

### `inputs/sample-sheet.csv`

One row per sequencing run. All samples here use the same guide library and genome.

| sample_id | expt_id | Run | genome | guides_filename | pam | scaffold | adapter_read1_5prime | adapter_read1_3prime |
|---|---|---|---|---|---|---|---|---|
| t0_rep1 | mtb_screen | SRRxxxxxxx | GCF_000195955.2 | guides.fasta | Sth1 | Sth1 | TGTAGCTTCTTTCGAGTACAAAAAC | CTCCCAGATTATATCTATCACTGATAGGGA |
| t0_rep2 | mtb_screen | SRRxxxxxxx | GCF_000195955.2 | guides.fasta | Sth1 | Sth1 | TGTAGCTTCTTTCGAGTACAAAAAC | CTCCCAGATTATATCTATCACTGATAGGGA |
| t6_rep1 | mtb_screen | SRRxxxxxxx | GCF_000195955.2 | guides.fasta | Sth1 | Sth1 | TGTAGCTTCTTTCGAGTACAAAAAC | CTCCCAGATTATATCTATCACTGATAGGGA |

`expt_id` groups samples for the fitness calculation — all samples to be compared must share the same value.

### `nextflow.config`

```nextflow
params {
    sample_sheet = "inputs/sample-sheet.csv"
    inputs       = "inputs"
    outputs      = "outputs"

    from_sra = true
    guides   = true
    rc       = true   // reverse complement guides before matching
}
```

### Run

```bash
nextflow run scbirlab/nf-crispriseq -r dev -c nextflow.config
```

### Output structure

```
outputs/
├── counts/                   # per-guide read counts, annotated with gene info
├── counts-by-experiment/     # per-experiment stacked tables (ready for fitness)
├── guides/                   # guide → genome GFF mappings
├── genome/                   # reference genome + eggNOG functional annotations
├── trimmed/                  # adapter-trimmed FASTQs
├── demultiplexed/            # guide-assigned FASTQs
├── plots/                    # count distribution plots
└── multiqc/                  # HTML QC report
```

The key file for most analyses is `outputs/counts/<sample_id>.annotated.tsv`: per-guide read counts joined to gene name, locus tag, and genomic coordinates.

---

## Processing steps

### Per genome (once per unique NCBI accession)

1. Download genome FASTA and GFF from NCBI.
2. **If `use_eggnog = true`:** annotate protein-coding genes with COG categories using [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper). Skipped by default as it requires downloading the eggNOG database (~GB-scale).
3. **If `guides = false`:** design all possible guide RNAs _de novo_ with [crispio](https://crispio.readthedocs.io) for the given PAM.
4. **If `guides = true`:** map provided guide sequences to the genome with `crispio map | crispio featurize` to assign each guide to a genomic feature.

### Per FASTQ file

1. **Trim** to adapter boundaries with `cutadapt` (quality + length filter).
2. **Extract UMIs** (if `use_umis = true`) with `umitools extract`. Optionally whitelist clone barcodes first (`use_clone_bc = true`).
3. **Count guides** — exact matching by default; set `allow_guide_errors = <n>` for mismatch-tolerant matching via `cutadapt`.
4. **Count reads** (and UMIs if applicable) per guide with `umitools count_tab`.
5. **Annotate** counts with the guide → gene table.
6. **Stack** per-sample tables into a per-experiment table.
7. **Plot** count distributions.

### QC

`fastqc` on raw reads; `multiqc` collates all logs into an HTML report.

### Fitness (optional)

If `do_fitness = true`, [`bartab fit`](https://github.com/scbirlab/bartab) models guide frequency change over the time course. See [Fitness calculation](#fitness-calculation).

---

## Requirements

You need Nextflow ≥23.10.0. The pipeline runs in containers by default — no other software installation required.

### On a local machine (Docker)

Install [Docker](https://docs.docker.com/get-docker/) and Nextflow:

```bash
conda install -c bioconda nextflow
# or: brew install nextflow
```

Docker is detected automatically. The container image is pulled on first run.

### On an HPC cluster (Singularity)

Most HPC systems provide Nextflow and Singularity as modules. Load them before running:

```bash
module load Singularity Nextflow
```

At the Crick, use the `standard` profile (see [Running on a cluster](#running-on-a-cluster)), which configures SLURM submission and Singularity automatically.

### Conda (alternative)

If you can't use containers, the pipeline falls back to conda. Set `conda.enabled = true` in your `nextflow.config` and ensure `conda` is available:

```nextflow
conda.enabled = true
```

### Setting `NXF_HOME`

If Nextflow can't write to its default cache location:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

---

## Inputs

### Parameters

**Required:**

| Parameter | Description |
|---|---|
| `sample_sheet` | Path to the CSV sample sheet |
| `inputs` | Directory containing guide files referenced in the sample sheet |
| `fastq_dir` | Path to local FASTQ directory. Required when `from_sra = false` (the default); omit if `from_sra = true`. |

**Optional (defaults shown):**

| Parameter | Default | Description |
|---|---|---|
| `outputs` | `"outputs"` | Output directory |
| `from_sra` | `false` | Pull FASTQs from SRA |
| `guides` | `false` | Provide a guide library; if `false`, guides are designed _de novo_ |
| `use_eggnog` | `false` | Run eggNOG-mapper to annotate genome with COG categories (requires large DB download) |
| `rc` | `false` | Reverse complement guide sequences before matching |
| `use_umis` | `false` | Reads contain UMIs |
| `use_clone_bc` | `false` | Reads contain clone barcodes (requires `use_umis = true`) |
| `allow_guide_errors` | `false` | Allow mismatches; set to an integer (e.g. `1`) to enable |
| `trim_qual` | `10` | Minimum Phred score for 3′ quality trimming |
| `min_length` | `15` | Minimum post-trimming read length |
| `max_length` | `false` | Maximum post-trimming read length (no cap by default) |
| `retain_5prime` | `false` | Retain rather than remove the 5′ adapter |
| `keep_missing_3prime` | `false` | Keep reads lacking the 3′ adapter |
| `keep_missing_5prime` | `false` | Keep reads lacking the 5′ adapter |
| `name_column` | `"Name"` | Guide CSV column containing guide names |
| `sequence_column` | `"guide_sequence"` | Guide CSV column containing sequences |
| `guide_length` | `20` | Spacer length for _de novo_ guide design |
| `do_fitness` | `false` | Run `bartab` fitness modelling |

Parameters can be set in `nextflow.config` or passed directly:

```bash
nextflow run scbirlab/nf-crispriseq -r dev \
    --sample_sheet /path/to/sample-sheet.csv \
    --inputs /path/to/inputs \
    --fastq_dir /path/to/fastqs \
    --guides --rc \
    --trim_qual 15 --min_length 90
```

---

### Sample sheet

A CSV with one row per sequencing sample.

#### Always required

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier |
| `expt_id` | Experiment identifier — samples sharing this value are grouped for fitness. **Required**: the pipeline will error on any row missing this value. |
| `genome` | [NCBI assembly accession](https://www.ncbi.nlm.nih.gov/datasets/genome/) (e.g. `GCF_000195955.2`) |
| `pam` | dCas9 PAM name (`Spy`, `Sth1`) or sequence (`NGG`, `NGRVAN`) |
| `scaffold` | sgRNA scaffold: `PerturbSeq` or `Sth1` |
| `adapter_read1_5prime` | 5′ adapter on R1 in [cutadapt format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequence left of the adapter is removed; the adapter itself is retained. |
| `adapter_read1_3prime` | 3′ adapter on R1. The adapter and everything to its right is removed. |

#### If using local FASTQs (`from_sra = false`, the default)

| Column | Description |
|---|---|
| `reads` | A substring of the FASTQ filename(s) in `fastq_dir`. The pipeline matches files as `*<reads>*`, so provide a unique fragment of the filename (e.g. a sample name or barcode), not a full glob. For paired-end data, the pattern should match **only R1**. |

#### If pulling from SRA (`from_sra = true`)

| Column | Description |
|---|---|
| `Run` | SRA Run accession |

#### If using paired-end reads

| Column | Description |
|---|---|
| `adapter_read2_5prime` | 5′ adapter on R2 |
| `adapter_read2_3prime` | 3′ adapter on R2 |

#### If using UMIs (`use_umis = true`)

| Column | Description |
|---|---|
| `umi_read1` | [umitools regex pattern](https://umi-tools.readthedocs.io/en/latest/regex.html) for UMI on R1 |
| `umi_read2` | umitools regex pattern for UMI on R2 (if applicable) |

#### If providing a guide library (`guides = true`)

| Column | Description |
|---|---|
| `guides_filename` | Filename of a guide CSV or FASTA inside `inputs` |

#### Example: SRA, single-end, no UMIs

| sample_id | expt_id | Run | genome | guides_filename | pam | scaffold | adapter_read1_5prime | adapter_read1_3prime |
|---|---|---|---|---|---|---|---|---|
| t0_rep1 | screen1 | SRRxxxxxxx | GCF_000195955.2 | guides.fasta | Sth1 | Sth1 | TGTAGCTTCTTTCGAGTACAAAAAC | CTCCCAGATTATATCTATCACTGATAGGGA |
| t6_rep1 | screen1 | SRRxxxxxxx | GCF_000195955.2 | guides.fasta | Sth1 | Sth1 | TGTAGCTTCTTTCGAGTACAAAAAC | CTCCCAGATTATATCTATCACTGATAGGGA |

#### Example: local FASTQs, paired-end, with UMIs

| sample_id | expt_id | reads | genome | guides_filename | pam | scaffold | adapter_read1_5prime | adapter_read1_3prime | umi_read1 |
|---|---|---|---|---|---|---|---|---|---|
| lib001 | screen2 | `FAU6865A42_*_R1` | GCA_003076915.1 | guides.csv | Spy | PerturbSeq | `^N{8}TCGACTGAGCTGAAAGAAT` | GTTTAAGAGCTATGCTGG | `^(?P<umi_1>.{8})(?P<discard_1>.{86}).+$` |

---

### Guide file formats

Three formats are accepted:

**CSV** (`.csv`):

| Name | guide_sequence |
|---|---|
| TargetA_guide1 | TCGACTGAGCTGAAAGAAT |
| TargetA_guide2 | GTTTAAGAGCTATGCTGGT |

**TSV** (`.tsv` or `.txt`) — same structure, tab-separated.

**FASTA** (any other extension):

```
>TargetA_guide1
TCGACTGAGCTGAAAGAAT
>TargetA_guide2
GTTTAAGAGCTATGCTGGT
```

Column names in the CSV/TSV can be changed with `name_column` and `sequence_column`. The pipeline maps guides to the genome and annotates targeted genes automatically — no need to provide coordinates.

### _De novo_ guide design

With `guides = false`, the pipeline uses [crispio](https://crispio.readthedocs.io) to design all possible guide RNAs for the given genome and PAM. Set `guide_length` (default 20 nt) to control spacer length.

---

## Outputs

| Path | Contents |
|---|---|
| `counts/<sample_id>.annotated.tsv` | Per-guide read counts joined to gene name, locus tag, coordinates |
| `counts-by-experiment/<expt_id>.counts.tsv.gz` | All samples in an experiment stacked into one table |
| `guides/` | Guide → genome GFF mappings |
| `genome/` | Reference genome FASTA/GFF; eggNOG annotations if `use_eggnog = true` |
| `trimmed/` | Trimmed FASTQs and cutadapt logs |
| `demultiplexed/` | Guide-assigned FASTQs |
| `plots/` | Count distribution histograms and correlations |
| `fitness/` | Fitness scores and plots (only with `do_fitness = true`) |
| `multiqc/` | HTML QC report |

---

## Fitness calculation

Set `do_fitness = true` to run [`bartab fit`](https://github.com/scbirlab/bartab) on the stacked count table for each experiment.

### Model types

Two models are available, selected automatically:

- **WLS** (default) — weighted least-squares fit to guide frequency over time.
- **HillFitnessModel** — dose-response model for concentration series; activated when `concentration_column` is set.

Normalisation uses either an external growth measurement supplied via `growth_column` (e.g. OD readings or generation counts from the sample sheet) or the read frequency of spike-in guides (`use_spike = true`). The `reference_guide` parameter specifies non-targeting control guides used to compute relative fitness after normalisation.

### Additional sample sheet columns for fitness

| Column | Description |
|---|---|
| `timepoint` | Numeric timepoint (generations, hours, etc.). Column name set by `timepoint_column`. |
| `culture_id` | Replicate culture identifier. Column name set by `culture_column`. |

### Fitness parameters

| Parameter | Default | Description |
|---|---|---|
| `reference_guide` | `"empty-vector"` | Name prefix of non-targeting control guides; used to compute relative fitness after normalisation |
| `timepoint_column` | `"timepoint"` | Sample sheet column with timepoint values |
| `culture_column` | `"culture_id"` | Sample sheet column identifying replicate cultures |
| `concentration_column` | `false` | Column with drug concentrations (activates HillFitnessModel) |
| `growth_column` | — | Sample sheet column containing growth measurements (OD, CFU, etc.) for normalisation |
| `growth` | `false` | Path to a separate TSV file mapping timepoints to growth measurements (alternative to `growth_column`) |
| `growth_type` | `"density"` | Units of the growth measurement: `"density"` (e.g. OD) or `"generations"` |
| `use_spike` | `false` | Normalise by spike-in guide frequency instead of a growth measurement. Only applied if `growth` is not set. |
| `negative` | `"ctrl_"` | Name prefix of negative/non-targeting control guides; shown distinctively in plots |
| `highlight_guides` | `false` | Comma-separated list of guide names to highlight in fitness plots |

### Example config with fitness

```nextflow
params {
    sample_sheet = "inputs/sample-sheet.csv"
    inputs       = "inputs"
    outputs      = "outputs"
    from_sra     = true
    guides       = true
    rc           = true

    do_fitness       = true
    reference_guide  = "empty-vector"
    timepoint_column = "generations"
    culture_column   = "replicate"
    growth_type      = "generations"
}
```

### Fitness outputs

| File | Description |
|---|---|
| `fitness/<expt_id>.bartab.h5ad` | Full results in AnnData format |
| `fitness/<expt_id>.bartab.csv` | Per-guide fitness parameters as a flat CSV |
| `fitness/plots/` | Fitness score plots |

---

## Running on a cluster

Load the required modules, then use the `standard` profile, which configures SLURM submission and Singularity automatically:

```bash
module load Singularity Nextflow
nextflow run scbirlab/nf-crispriseq -r dev -profile standard -c nextflow.config
```

The `standard` profile also enables email notification on completion (to `$USER@crick.ac.uk`) and writes a DAG of the pipeline graph.

Default resource allocations by process label:

| Label | CPUs | Memory | Time |
|---|---|---|---|
| (default) | 1 | 8 GB | 12 h |
| `big_cpu` | 16 | 16 GB | 12 h |
| `big_time` | 1 | 16 GB | 4 days |
| `some_mem` | 1 | 16 GB | 12 h |
| `med_mem` | 1 | 64 GB | 12 h |
| `big_mem` | 16 | 128 GB | 12 h |
| `gpu` | 2 | 128 GB | 4 h |
| `gpu_single` | 2 | 128 GB | 7 days |

Override in your `nextflow.config`:

```nextflow
process {
    withLabel: big_time {
        time   = '7d'
        memory = 32.GB
    }
}
```

To run locally (e.g. for testing with Docker):

```bash
nextflow run scbirlab/nf-crispriseq -r dev -profile local -c nextflow.config
```

---

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-crispriseq/issues).

---

## Further help

- [bartab](https://github.com/scbirlab/bartab) — fitness modelling
- [crispio](https://crispio.readthedocs.io/en/stable/index.html) — guide design and genome mapping
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) — adapter trimming
- [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) — functional annotation
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
- [umi-tools](https://umi-tools.readthedocs.io/en/latest/index.html)
