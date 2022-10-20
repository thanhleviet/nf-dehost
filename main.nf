#!/usr/bin/env nextflow


process FASTQC {
    
    label "CHANGE_ME"
    
    tag {sample_id}
    
    cpus 4

    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("fastqc_${sample_id}_logs"), emit: logs
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs $reads -q ${reads}
    """
}


process FASTP {
    
    label "CHANGE_ME"
    
    tag {sample_id}
    
    cpus 16

    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(simple_sample_id), path("*.json"), emit: logs
    script:
    simple_sample_id = file(reads[0]).getBaseName() - ".fastq"

    """
    fastp -i ${reads[0]} \
    -o /dev/null \
    -I ${reads[1]} \
    -O 1>/dev/null \
    -j \
    -w ${task.cpus}
    mv fastp.json ${simple_sample_id}.fastp.json
    """
}

process BOWTIE2 {
    
    label "CHANGE_ME"
    
    tag {sample_id}
    
    cpus 8

    input:
    tuple val(sample_id), path(reads)
    path(ref)

    output:
    tuple val(sample_id), path("*.non_host.fastq.gz"), emit: non_host

    script:
    ref_basename = ref[0].getSimpleName()
    """
    bowtie2 --very-sensitive-local -x ${ref_basename} -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} | samtools fastq -@ ${task.cpus} -1 ${sample_id}_1.non_host.fastq -2 ${sample_id}_2.non_host.fastq -0 /dev/null -s /dev/null -n -f 4 -
    pigz -p ${task.cpus} *.non_host.fastq
    """
}

process MULTIQC {
    
    label "CHANGE_ME"
    
    tag {"Running"}
    
    cpus 4

    input:
    path(fastqc_logs)

    output:
    path("multiqc_*"), emit: report

    script:
    """
    multiqc .
    """
}

workflow FASTP_PRE_DEHOST {
    take:
    reads
    main:
    FASTP(reads)
    emit:
    logs = FASTP.out.logs
}

workflow FASTP_POST_DEHOST {
    take:
    reads
    main:
    FASTP(reads)
    emit:
    logs = FASTP.out.logs
}

ch_reads = Channel.fromPath(params.sample_sheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, [row.sr1, row.sr2]) }

ch_ref = Channel.fromPath(params.ref + "*.bt2", checkIfExists: true)
                .collect()

log.info("HUMAN READS REMOVAL")
log.info("Sample sheet: " + params.sample_sheet)
log.info("Database    : " + params.ref)
log.info("Outdir      : " + params.outdir)

workflow {
    FASTP_PRE_DEHOST(ch_reads)
    BOWTIE2(ch_reads, ch_ref)
    FASTP_POST_DEHOST(BOWTIE2.out.non_host)
    ch_multiqc = FASTP_PRE_DEHOST.out.logs.map{it[1]}.collect().combine(FASTP_POST_DEHOST.out.logs.map{it[1]}.collect())
    MULTIQC(ch_multiqc)
}
