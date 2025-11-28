process bam_to_fastq {
    memory = { 8.GB * task.attempt }
    cpus 2
    errorStrategy = { (task.exitStatus == 143 || task.exitStatus == 137) ? 'retry' : 'terminate' }
    maxRetries 2
    
    //container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'
    
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.fastq")
    
    script:
    """
    samtools fastq -T '*' ${bam} > ${sam}_${type}.fastq
    """
}

process minimap2_ont {
    memory = { 32.GB * task.attempt }
    cpus 8
    errorStrategy = { (task.exitStatus == 143 || task.exitStatus == 137) ? 'retry' : 'terminate' }
    maxRetries 2
    
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(fastq)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam")
    
    script:
    """
    minimap2 -x map-ont -a -t ${task.cpus} --MD ${params.ref_fasta} ${fastq} | \\
      samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    """
}

process minimap2_pacbio {
    memory = { 32.GB * task.attempt }
    cpus 8
    errorStrategy = { (task.exitStatus == 143 || task.exitStatus == 137) ? 'retry' : 'terminate' }
    maxRetries 2
    
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(fastq)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam")
    
    script:
    """
    minimap2 -x map-hifi -a -t ${task.cpus} --MD ${params.ref_fasta} ${fastq} | \\
      samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    """
}

process index_bam {
    memory = { 4.GB * task.attempt }
    cpus 4
    errorStrategy = { (task.exitStatus == 143 || task.exitStatus == 137) ? 'retry' : 'terminate' }
    maxRetries 2
   
    tag { sam }
    publishDir "output/aligned_bams", mode: 'copy'
    
    input:
    tuple val(sam), val(type), path(bam)
    
    output:
    tuple val(sam), val(type), path(bam), path("${bam}.bai")
    
    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

// Example workflow showing how to connect the processes
workflow minimap2_ubam_ont {
    take:
    input_bams  // channel: tuple val(sam), val(type), path(bam)
    
    main:
    // Step 1: Convert BAM to FASTQ
    bam_to_fastq(input_bams)
    
    // Step 2: Align with minimap2 and sort
    minimap2_ont(
        bam_to_fastq.out)
    
    // Step 3: Index the sorted BAM
    index_bam(minimap2_ont.out)
    
    emit:
    aligned_bams = index_bam.out  // tuple val(sam), val(type), path(bam), path(bai)
}

workflow minimap2_ubam_pacbio {
    take:
    input_bams  // channel: tuple val(sam), val(type), path(bam)
    
    main:
    // Step 1: Convert BAM to FASTQ
    bam_to_fastq(input_bams)
    
    // Step 2: Align with minimap2 and sort
    minimap2_pacbio(
        bam_to_fastq.out)
    
    // Step 3: Index the sorted BAM
    index_bam(minimap2_pacbio.out)
    
    emit:
    aligned_bams = index_bam.out  // tuple val(sam), val(type), path(bam), path(bai)
}


process minimap2_ubam_illumina {
    memory = { 32.GB * task.attempt }
    errorStrategy = { (task.exitStatus == 143 || task.exitStatus == 137) ? 'retry' : 'terminate' }

    tag {sam}
    cpus 8
    input:
    tuple val(sam), val(type), path(bam)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai")
    
    script:
    """
    samtools fastq -T "*" -@ ${task.cpus} ${bam} | \\
        minimap2 -ax sr -t ${task.cpus} ${params.ref_fasta} - | \\
        samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    samtools index -@ ${task.cpus} ${sam}_${type}.sorted.bam
    """
}
