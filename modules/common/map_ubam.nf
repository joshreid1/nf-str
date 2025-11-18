process minimap2_ubam_ont {
    memory = { 32.GB * task.attempt }
    errorStrategy = { (task.exitStatus == 143 || task.exitStatus == 137) ? 'retry' : 'terminate' }

    cpus 8
    tag { sam }
    

    input:
    tuple val(sam), val(type), path(bam)

    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai")
    
    script:
    """
    samtools fastq -T '*' ${bam} | \\
      minimap2 -x map-ont -a  -t ${task.cpus} --MD  ${params.ref_fasta} - | \\
      samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    samtools index -@ ${task.cpus} ${sam}_${type}.sorted.bam
    """
}

process minimap2_ubam_pacbio {
    memory = { 32.GB * task.attempt }
    errorStrategy = { (task.exitStatus == 143 || task.exitStatus == 137) ? 'retry' : 'terminate' }

    cpus 8
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)

    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai")
    
    script:
    """
    samtools fastq -T '*' ${bam} | \\
      minimap2 -x map-hifi -a -t ${task.cpus} --MD  ${params.ref_fasta} - | \\
      samtools sort -O bam -@ ${task.cpus} -o ${sam}_${type}.sorted.bam -
    samtools index -@ ${task.cpus} ${sam}_${type}.sorted.bam
    """
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
