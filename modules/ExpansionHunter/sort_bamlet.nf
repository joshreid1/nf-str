process sort_bamlet {
    cpus 1
    memory '1 GB'
    time '1 hours'
    publishDir "output/eh5/", mode: "copy"
    tag { sam }
    
    input:
    tuple val(sam), path(bamlet), path(json), path(vcf)

    output:
    tuple val(sam), path("${sam}_sorted_realigned.bam"), path("${sam}_sorted_realigned.bam.bai"), path(vcf)
    
    script:
    """
    samtools sort -o ${sam}_sorted_realigned.bam ${bamlet}
    samtools index ${sam}_sorted_realigned.bam
    """
}