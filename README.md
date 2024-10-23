# 🧮 CRISPRi-seq pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-crispriseq/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**scbirlab/nf-crispriseq** is a Nextflow pipeline to count and annotate guide RNAs in demultiplexed 
FASTQ files, optionally with UMIs, and optionally modelling fitness changes.

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

Per genome:

0. If no guide RNAs are provided, generate all possible guide RNAs from the provided genome
1. For each input set of guide RNAs, map to the reference genome and GFF to get genomic feature annotations

Per FASTQ file:

1. Filter and trim reads to adapters using `cutadapt`. This ensures reads used downstream 
have the expected features and are trimmed so that the features are in predictable places.
2. (Optionally) Extract UMIs using `umitools extract`.
3. Find guide RNA matches using `cutadapt`.
4. Count UMIs (if using) and reads per guide RNA using `umitools count_tab`.
7. Plot histograms and correlations of UMI and read count distributions.

Optionally [work in progress]:

8. If the data are from a time-course, calculate fitness per guide RNA, per condition.

### Other steps

1. Get FASTQ quality metrics with `fastqc`.
2. Compile the logs of processing steps into an HTML report with `multiqc`.

## Requirements

### Software

You need to have Nextflow and `conda` installed on your system.

#### First time using Nextflow?

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow
```

Otherwise, if it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet)  and, optionally, a [`nextflow.config` file](#inputs) 
in the directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-crispriseq
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which 
will not be automatically updated. If you want to ensure that you're running a version of 
the pipeline, use the `-r <version>` flag. For example,

```bash 
nextflow run scbirlab/nf-crispriseq -r v0.0.1
```

A list of versions can be found by running `nextflow info scbirlab/nf-crispriseq`.

For help, use `nextflow run scbirlab/nf-crispriseq --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` 
will be installed. This can take around 10 minutes.

If your run is unexpectedly interrupted, you can restart from the last completed step using the `-resume` flag.

```bash 
nextflow run scbirlab/nf-crispriseq -resume
```

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV containing sample IDs matched with FASTQ filenames, references, 
and adapter sequences
- `fastq_dir`: path to directory containing the FASTQ files (optionally GZIPped)
- `inputs`: path to directory containing files referenced in the `sample_sheet`, such as lists of guide RNAs.

The following parameters have default values that can be overridden if necessary.

- `output = "outputs"`: path to directory to put output files
- `sample_names = "sample_id"`: column of `sample_sheet` to take as a sample identifier. Use `"Run"` for SRA table inputs.
- `use_umis = false`: Whether to the reads include UMIs
- `from_sra = false`: Whether the FASTQ files should be pulled from the SRA instead of provided as local files
- `guides = true`: Whether the name of a CSV of guide sequences is provided in the `sample_sheet`
- `name_column = "Name"`: If using a CSV of guide RNA sequences (`guides = true`), the column containing the name of each guide
- `sequence_column = "guide_sequence"`: If using a CSV of guide RNA sequences (`guides = true`), the column containing the sequence of each guide
- `rc = false`: Whether to reverse complement the guide sequences before mapping.
- `trim_qual = 10` : For `cutadapt`, the minimum Phred score for trimming 3' calls
- `min_length = 105` : For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
   
    sample_sheet = "/path/to/sample-sheet.csv"
    fastq_path = "/path/to/fastqs"
    guides = "/path/to/reference"

    // Optional
    rc = true
    guides = true
    trim_qual = 15
    min_length = 90

}
```

Alternatively, you can provide these on the command line:

```bash
nextflow run scbirlab/nf-crispriseq -r v0.0.1 \
    --sample_sheet /path/to/sample_sheet.csv \
    --fastq /path/to/fastqs \
    --reference /path/to/reference \
    --rc --guides \
    --trim_qual 15 --min_length 90
``` 

### Sample sheet

The sample sheet is a CSV file providing information about each sample: which FASTQ files belong 
to it, the reference genome accession number, adapters to be trimmed off, (optionally) the UMI  
scheme, (optionally) the name of a table of known guide RNAs, and (optionally) experimental conditions if calculating fitness.

The file must have a header with the column names below, and one line per sample to be processed.

- `sample_id`: the unique name of the sample
- `genome`: The [NCBI assembly accession](https://www.ncbi.nlm.nih.gov/datasets/genome/) number for the organism that the guide RNAs are targeting. This number starts with "GCF_" or "GCA_".
- `pam`: The name (e.g. "Spy" or "Sth1") or sequence (e.g "NGG" or "NGRVAN") of the dCas9 PAM used in the experiment
- `scaffold`: The name of the sgRNA scaffold ("PerturbSeq" or "Sth1") used in the experiment.
The pipeline will look for files matching `<fastq_dir>/*<dir>*`, and should match **only the forward read** if you had paired-end sequencing.
- `adapter_5prime`: the 5' adapter on the forward read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequence _to the left_ will be removed, but the adapters themselves will be retained.
- `adapter_3prime`: the 3' adapter on the forward read to trim to in [`cutadapt` format](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences). Sequence _**matching the adapter** and everything to the right_ will be removed.

If you have set `from_sra = false` (the default):
- `reads`: the search glob to find FASTQ files for each sample in `fastq_dir` (see [config](#inputs)). 
Otherwise with `from_sra = true`:
- `Run`: the SRA Run ID

If you have set `use_umis = true` (the default):
- `umi_pattern`: the cell barcode and UMI pattern in [`umitools` regex format](https://umi-tools.readthedocs.io/en/latest/regex.html#regex-regular-expression-mode) for the forward read

If you have set `guides = true` (the default):
- `guides_filename`: the name of a file in the inputs directory containing guide sequences. 

Here is an example of the sample sheet:

| sample_id | genome          | reads           | guides_filename | pam  | scaffold   | adapter_5prime           | adapter_3prime     | umi_pattern                             | 
| --------- | --------------- | --------------- | --------------- | ---- | ---------- | ------------------------ | ------------------ | --------------------------------------- |  
| lib001    | GCA_003076915.1 | FAU6865A42_*_R1 | guides.csv      | Spy  | PerturbSeq | ^N{8}TCGACTGAGCTGAAAGAAT | GTTTAAGAGCTATGCTGG | ^(?P<umi_1>.{8})(?P<discard_1>.{86}).+$ | 
| lib002    | GCA_003076915.1 | FAU6865A43_*_R1 | guides.csv      | Spy  | PerturbSeq | ^N{8}TCGACTGAGCTGAAAGAAT | GTTTAAGAGCTATGCTGG | ^(?P<umi_1>.{8})(?P<discard_1>.{86}).+$ | 

And here is an example of the `guides_filename` (`guides.csv` in this example):

| Name       | guide_sequence      |
| ---------- | ------------------- |
| guide001   | TCGACTGAGCTGAAAGAAT |
| guide002   | GTTTAAGAGCTATGCTGGT |

It is also possible to provide the guides as a fasta file:

```
>guide001
TCGACTGAGCTGAAAGAAT
>guide002
GTTTAAGAGCTATGCTGGT
```

You don't need to provide gene anntotation infomation, because the pipeline will map these guides back to the genome and annotate the features for you.

## Outputs

Outputs are saved in the `output` defined in the [config file](#inputs). They are organised under 
three directories:

- `processed`: FASTQ files and logs resulting from trimming and UMI extraction
- `mapped`: FASTQ files and logs resulting mapping features
- `counts`: tables and plots relating to UMI and read counts
- `multiqc`: HTML report on processing steps

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-crispriseq/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [crispio](https://crispio.readthedocs.io/en/stable/index.html)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
- [python](https://www.python.org/doc/)
- [matplotlib](https://matplotlib.org/stable/)
- [umi-tools](https://umi-tools.readthedocs.io/en/latest/index.html)