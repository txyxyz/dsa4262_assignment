#!/usr/bin/env nextflow

// Adapted from the SG-NEx long-read RNA-seq workshop workflow:
// GoekeLab/sg-nex-data/docs/colab/workflow_longReadRNASeq.nf

params.refFa = '/path/to/ref.fa'
params.refGtf = '/path/to/ref.gtf'
params.reads = '/path/to/*.fastq.gz'
params.outdir = 'results'
params.use_annotation = true


// ------------------------------------------------------------
// 1. Alignment with Minimap2
// ------------------------------------------------------------

process MINIMAP2_ALIGN {

    tag "${sample_id}"
    maxForks 1

    input:
    tuple val(sample_id), path(reads)
    path refFa

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    def minimap_opts = sample_id.contains("directRNA") ?
        "-ax splice -uf -k14" :
        "-ax splice"

    """
    minimap2 ${minimap_opts} ${refFa} ${reads} > ${sample_id}.sam
    """
}


// ------------------------------------------------------------
// 2. SAM to BAM conversion
// ------------------------------------------------------------

process SAM_TO_BAM {

    tag "${sample_id}"

    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_sam)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    """
    samtools view -b ${reads_sam} > ${sample_id}.bam
    """
}


// ------------------------------------------------------------
// 3. Quality control
// ------------------------------------------------------------

process QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_bam)

    output:
    path "${sample_id}.flagstat.txt"

    script:
    """
    samtools flagstat ${reads_bam} > ${sample_id}.flagstat.txt
    """
}


// ------------------------------------------------------------
// 4. Transcript discovery and quantification with Bambu
// ------------------------------------------------------------

process BAMBU {

    publishDir "${params.outdir}/bambu", mode: 'copy'

    input:
    path refFa
    path refGtf
    path reads_bam

    output:
    path "counts_transcript.txt"
    path "counts_gene.txt"
    path "extended_annotations.gtf"

    script:
    def annotation_code = params.use_annotation ?
        """
        annotations <- prepareAnnotations("${refGtf}")
        se <- bambu(
            reads = c(${reads_bam.collect { "\"${it}\"" }.join(", ")}),
            annotations = annotations,
            genome = "${refFa}",
            ncore = 2
        )
        """ :
        """
        se <- bambu(
            reads = c(${reads_bam.collect { "\"${it}\"" }.join(", ")}),
            genome = "${refFa}",
            NDR = 1,
            ncore = 2
        )
        """

    """
    #!/usr/bin/env Rscript --vanilla

    library(bambu)

    ${annotation_code}

    writeBambuOutput(se, path = "./")
    """
}


// ------------------------------------------------------------
// Workflow
// ------------------------------------------------------------

workflow {

    reads_ch = Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { read ->
	    def sample_id = read.name.replaceFirst(/\.fastq\.gz$/,'')
            tuple(sample_id, read)
        }

    ref_ch = Channel.value(file(params.refFa))
    gtf_ch = Channel.value(file(params.refGtf))

    MINIMAP2_ALIGN(reads_ch, ref_ch)

    SAM_TO_BAM(MINIMAP2_ALIGN.out)

    // Split BAM channel so both QC and Bambu can consume it
    bam_for_qc = SAM_TO_BAM.out.bam.map { sample_id, bam ->
        tuple(sample_id, bam)
    }

    bam_for_bambu = SAM_TO_BAM.out.bam
        .map { sample_id, bam -> bam }
        .collect()

    QC(bam_for_qc)

    BAMBU(ref_ch, gtf_ch, bam_for_bambu)
}
