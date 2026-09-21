# Long-read RNA-seq Nextflow Workflow

This repository contains the Nextflow workflow used for Task 5 of the genomics assignment.

The workflow was adapted from the long-read RNA-seq workflow provided in the Introduction to Genomics 3 workshop, originally from the SG-NEx project.

## Workflow

The pipeline performs the following steps:

1. Alignment of long-read RNA-seq reads to the reference genome using Minimap2.
2. Conversion of SAM files to BAM files using Samtools.
3. Quality control using `samtools flagstat`.
4. Transcript discovery and gene/transcript quantification using Bambu.

The workflow processes both direct RNA and cDNA sequencing data. Protocol-specific Minimap2 parameters are used:

- Direct RNA: `-ax splice -uf -k14`
- cDNA: `-ax splice`

## Input

The workflow requires:

- FASTQ files containing long-read RNA-seq reads
- Reference genome in FASTA format
- Genome annotation in GTF format

## Output

The workflow produces:

- BAM alignment files
- `samtools flagstat` QC reports
- `extended_annotations.gtf`
- `counts_transcript.txt`
- `counts_gene.txt`

## Running the workflow

### Scenario 1: Bambu with genome annotations

```bash
nextflow run workflow_longReadRNASeq_task5.nf \
  --reads "/path/to/*.fastq.gz" \
  --refFa "/path/to/reference.fa" \
  --refGtf "/path/to/annotation.gtf" \
  --use_annotation true \
  --outdir "results_with_annotation" \
  -with-report report_with_annotation.html
````

### Scenario 2: Bambu without genome annotations
```bash
nextflow run workflow_longReadRNASeq_task5.nf \
  --reads "/path/to/*.fastq.gz" \
  --refFa "/path/to/reference.fa" \
  --refGtf "/path/to/annotation.gtf" \
  --use_annotation false \
  --outdir "results_without_annotation" \
  -resume \
  -with-report report_without_annotation.html
````

## Software
- Nextflow 26.04.6
- Bambu 3.12.1
- Minimap2
- Samtools
- R 4.5.3

## Acknowledgement

This workflow was adapted from the long-read RNA-seq Nextflow workflow provided in the Introduction to Genomics 3 workshop and the SG-NEx project.

Original SG-NEx repository:

https://github.com/GoekeLab/sg-nex-data
