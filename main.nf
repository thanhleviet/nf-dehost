#!/usr/bin/env nextflow


process FASTQC {
    
    label "CHANGE_ME"
    
    conda "bioconda::fastqc"
    
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
    
    conda "bioconda::fastp"

    tag {sample_id}
    
    cpus 5

    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(simple_sample_id), path("*.json"), emit: logs
    script:
    simple_sample_id = file(reads[0]).getBaseName() - ".fastq"

    """
    fastp -i ${reads[0]} \
    -I ${reads[1]} \
    -j fastp.json \
    -w ${task.cpus}
    mv fastp.json ${simple_sample_id}.fastp.json
    """
}

process BOWTIE2 {
    
    label "CHANGE_ME"
    conda "bioconda::bowtie2=2.5.1 bioconda::samtools"
    
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
    
    conda "bioconda::multiqc=1.14.0"
    
    tag {"Running"}
    
    cpus 4

    input:
    path(fastqc_logs)

    output:
    path("multiqc_*"), emit: report
    path("sample_reads.tsv")

    shell:
    """
    for f in *.json; do echo \$f,\$(jq .summary.after_filtering.total_reads \$f) | sed 's/.fastp.json//g';done > sample_reads.tsv
    multiqc .
    """
}

process fastqMergeLanes {
  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: 'copy'

  tag { sampleName }
  cpus 8 
  
  input:
  tuple val(sampleName), path(forward), path(reverse)
  
  output:
  tuple val(sampleName), path("${sampleName}_R1.fastq.gz"), path("${sampleName}_R2.fastq.gz")
  
  script:
  """
  zcat $forward | pigz -p ${task.cpus} - > ${sampleName}_R1.fastq.gz
  zcat $reverse | pigz -p ${task.cpus} - > ${sampleName}_R2.fastq.gz
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

workflow wf_fastqMergeLanes {
  take:
    ch_filePairs
  
  main:
    // Group 4 lane elements into 1 element in the channel
    ch_filePairs
            .map {
            it -> [it[0].replaceAll(~/\_L00[1,2,3,4]/,""), it[1], it[2]]
            }
            .groupTuple(by:0)
            .set { ch_reads_four_lanes }

    fastqMergeLanes(ch_reads_four_lanes)
  
  emit: fastqMergeLanes.out
}

ch_input = Channel.fromPath(params.sample_sheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, [row.R1, row.R2]) }

ch_ref = Channel.fromPath(params.ref + "*.bt2", checkIfExists: true)
                .collect()


log.info("HUMAN READS REMOVAL")
log.info("Sample sheet: " + params.sample_sheet)
log.info("Database    : " + params.ref)
log.info("Outdir      : " + params.outdir)

workflow {
    if (params.merge_lanes) {
        wf_fastqMergeLanes(ch_input)
        ch_reads = wf_fastqMergeLanes.out.map {it -> tuple(it[0],[it[1],it[2]])}
    } else {
        ch_reads = ch_input
    }

    FASTP_PRE_DEHOST(ch_reads)
    BOWTIE2(ch_reads, ch_ref)
    FASTP_POST_DEHOST(BOWTIE2.out.non_host)
    ch_multiqc = FASTP_PRE_DEHOST.out.logs.map{it[1]}.collect().combine(FASTP_POST_DEHOST.out.logs.map{it[1]}.collect())
    MULTIQC(ch_multiqc)
}
