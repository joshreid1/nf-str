process sort_bamlet {
    cpus 1
    memory '1 GB'
    time '1 hours'
    publishDir "progress/eh5/", mode: "symlink"
    
    input:
    tuple val(sam), path(bamlet)
    
    output:
    tuple val(sam), path("${sam}_sorted_realigned.bam"), path("${sam}_sorted_realigned.bam.bai")
    
    script:
    """
    micromamba run -n base samtools sort -o ${sam}_sorted_relaigned.bam ${bamlet}
    micromamba run -n base tabix -p bam ${sam}_sorted_relaigned.bam
    """
}