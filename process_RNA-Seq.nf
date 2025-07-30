#!/usr/bin/env nextflow

// Create a channel of input FASTQ file pairs by matching R1 files and deriving R2 and sample name
Channel
    .fromPath(params.input)
    .map { read1 ->
        def read2 = file(read1.toString().replaceFirst('_1\\.fastq\\.gz$', '_2.fastq.gz'))
        def sample = read1.getBaseName().split('_')[0]
        tuple(sample, read1, read2)
    }
    .set { read_pairs }

// FASTP - Trims and filters input reads
process Fastp {
    tag "$sample"
    cpus 6

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_fastp_R1.fastq.gz"),
          path("${sample}_fastp_R2.fastq.gz"),
          path("${sample}_fastp_S1.fastq.gz"),
          path("${sample}_fastp_S2.fastq.gz"),
          path("${sample}_fastp.json"),
          path("${sample}_fastp.html"),
          path("${sample}_fastp.log")

   publishDir 'results/fastp_cleaned_reads', mode: 'copy', pattern: "*fastp_{R,S}*.fastq.gz"
   publishDir 'results/logs/fastp', mode: 'copy', pattern: "*.{json,html,log}"

    script:
    """
    fastp \
        -i $read1 -I $read2 \
        --out1 ${sample}_fastp_R1.fastq.gz \
        --out2 ${sample}_fastp_R2.fastq.gz \
        --unpaired1 ${sample}_fastp_S1.fastq.gz \
        --unpaired2 ${sample}_fastp_S2.fastq.gz \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_front_window_size 5 \
        --cut_front_mean_quality 20 \
        --correction \
        -l 30 \
        --thread ${task.cpus} \
        -j ${sample}_fastp.json \
        -h ${sample}_fastp.html \
        2> ${sample}_fastp.log
    """
}

// AlignReads - Maps reads to reference genome using HISAT2 then sorts
process AlignReads {
    tag "$sample"
    cpus 6

    input:
    tuple val(sample),
          path(r1), path(r2), path(s1), path(s2),
          path(json), path(html), path(log)

    output:
    tuple val(sample),
          path("${sample}_Halign_sort.bam"),
          path("${sample}_Halign_sort.bam.bai"),
          path("${sample}_align.summ.txt"),
          path("${sample}.unmapped.fq")

    publishDir 'results/bams', mode: 'copy', pattern: "*_Halign_sort.bam{,.bai}"
    publishDir 'results/logs/align', mode: 'copy', pattern: "*_align.summ.txt|*.unmapped.fq"

    script:
    """
    # Extract RGID (flowcell ID) and PU (lane) from first FASTQ read
    RG_LINE=\$(zcat $r1 | head -n 1)
    RGID=\$(echo \$RG_LINE | awk -F ':' '{print \$3}')
    PU=\$(echo \$RG_LINE | awk -F ':' '{print \$4}')

    echo "Extracted RGID: \$RGID"
    echo "Extracted PU: \$PU"
    hisat2 -p ${task.cpus} --mp 2,0 -q \
        -x ${params.ref} \
        -1 $r1 -2 $r2 \
        -U $s1,$s2 \
        --summary-file ${sample}_align.summ.txt \
        --un ${sample}.unmapped.fq \
    | samtools view -Sb - \
    | samtools sort -o ${sample}_Halign_sort.bam

    picard -Xmx20g AddOrReplaceReadGroups \
        I=${sample}_Halign_sort.bam \
        O=${sample}_Halign_sort.bam.tmp \
        SO=coordinate \
        RGID=\$RGID \
        LB=evans_MSU \
        PL=illumina \
        PU=\$PU \
        SM=${sample} \
        TMP_DIR=${params.tmp}

    mv ${sample}_Halign_sort.bam.tmp ${sample}_Halign_sort.bam
    samtools index ${sample}_Halign_sort.bam
    """
}

// QualimapQC - Performs quality control on BAM files using Qualimap
process QualimapQC {
    tag "$sample"
    cpus 6
    memory '25 GB'

    input:
    tuple val(sample), path(bam), path(bai), path(summary), path(unmapped)

    output:
    path("${sample}_bamqc")

    publishDir 'results/qualimap', mode: 'copy'

    script:
    """
    mkdir ${sample}_bamqc
    qualimap bamqc -bam $bam -outdir ${sample}_bamqc -outformat HTML -nt ${task.cpus}
    """
}

// MultiQC_Fastp - Aggregates Fastp reports using MultiQC
process MultiQC_Fastp {
    tag "multiqc_fastp"
    cpus 1

    input:
    path fastp_logs_files

    output:
    path "multiqc_fastp_report.html"

    publishDir "results/multiqc/fastp", mode: 'copy'

    script:
    """
    multiqc $fastp_logs_files --filename multiqc_fastp_report.html --outdir .
    """
}

// MultiQC_Qualimap - Aggregates Qualimap results using MultiQC
process MultiQC_Qualimap {
    tag "multiqc_qualimap"
    cpus 1

    input:
    path qualimap_dirs

    output:
    path "multiqc_qualimap_report.html"

    publishDir "results/multiqc/qualimap", mode: 'copy'

    script:
    """
    multiqc $qualimap_dirs --filename multiqc_qualimap_report.html --outdir .
    """
}

// FeatureCounts - Counts reads for all BAMs using featureCounts
process FeatureCounts {
    tag "featurecounts_combined"
    cpus 31
    memory '8 GB'

    input:
    path bam_files
    path gtf_file

    output:
    path "featureCounts_combined.txt"
    path "featureCounts_combined_log.txt"

    publishDir "results/featurecounts", mode: 'copy'

    script:
    """
    featureCounts -p --countReadPairs -O -F GTF \
        -T ${task.cpus} \
        -a $gtf_file \
        -o featureCounts_combined.txt \
        ${bam_files.join(' ')} 2> featureCounts_combined_log.txt
    """
}

// WORKFLOW block - Orchestrates the execution of processes
workflow {
    fastp_results = read_pairs | Fastp
    aligned = fastp_results | AlignReads
    qualimap_results = aligned | QualimapQC

    bam_files_ch = aligned.map { it[1] }.collect()
    gtf_ch = Channel.fromPath(params.gtf)
    FeatureCounts(bam_files_ch, gtf_ch)

    fastp_logs = fastp_results.flatMap { sample, r1, r2, s1, s2, json, html, log -> [json, html, log] }.collect()
    MultiQC_Fastp(fastp_logs)

    qualimap_dirs = qualimap_results.collect()
    MultiQC_Qualimap(qualimap_dirs)
}

// Cleanup work directory on successful completion
workflow.onComplete {
    if (workflow.success) {
        println "Workflow succeeded. Cleaning up work directory..."
        def dir = file(workDir)
        if (dir.exists()) {
            dir.listFiles().each { it.deleteDir() }
            println "Deleted contents of work directory: ${dir}"
        } else {
            println "Work directory not found: ${dir}"
        }
    } else {
        println "Workflow failed. Work directory preserved for debugging."
    }
}

