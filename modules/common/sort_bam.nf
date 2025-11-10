process sort_index_bam {
    cpus 8
    memory '10 GB'
    time '1 hours'
    tag { sam }
    
    input:
    tuple val(sam), path(bam)

    output:
    tuple val(sam), path("${sam}_sorted.bam"), path("${sam}_sorted.bam.bai")
    
    script:
    """
    samtools sort --threads 8 -o ${sam}_sorted.bam ${bam}
    samtools index --threads 8 ${sam}_sorted.bam
    """
}