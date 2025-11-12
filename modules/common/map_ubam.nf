process minimap2_ubam_ont {
    cpus 8
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)

    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai")
    
    script:
    """
    samtools fastq -T '*' ${bam} | \\
      minimap2 -x map-ont -a  -t ${task.cpus}--MD  ${params.ref_fa} - | \\
      samtools sort  -@ ${task.cpus} -o ${sam}_${type}.sorted.bam --write-index -
    """
}

process minimap2_ubam_pacbio {
    cpus 8
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)

    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai")
    
    script:
    """
    samtools fastq -T '*' ${bam} |  minimap2 -x map-hifi -a -t ${task.cpus} --MD  ${params.ref_fa} - | samtools sort -@ ${task.cpus} -o ${sam}_${type}.sorted.bam --write-index -
    """
}

process minimap2_ubam_illumina {
    tag { sam}
    cpus 8
    input:
    tuple val(sam), val(type), path(bam)
    
    output:
    tuple val(sam), val(type), path("${sam}_${type}.sorted.bam"), path("${sam}_${type}.sorted.bam.bai")
    
    script:
    """
    samtools fastq -@ ${task.cpus} ${bam} | \\
        minimap2 -ax sr -t ${task.cpus} ${params.ref_fa} - | \\
        samtools sort -@ ${task.cpus} -o ${sam}.sorted.bam --write-index -

    """
}