process index_bam {
    cpus 8
    memory '10 GB'
    time '1 hours'
    tag { sam }
    
    input:
    tuple val(sam), val(type), path(bam)

    output:
    tuple val(sam), val(type),  path(bam), path("${bam}*")
    
    script:
    """
    samtools index --threads 8 ${bam} 
    """
}

process sort_bam {
    cpus 8
    memory '10 GB'
    time '1 hours'
    publishDir "output/eh5/", mode: "copy"
    tag { sam }
    
    input:
    tuple val(sam), path(bam)

    output:
    tuple val(sam), path("*_sorted.bam")
    
    script:
    def bam_sorted = bam.replaceAll('.bam', '_sorted.bam')
    """
    samtools sort --threads 8 -o ${bam_sorted} ${bam}
    """
}
