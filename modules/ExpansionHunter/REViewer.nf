params.catalog = "${projectDir}/catalogues/ExpansionHunter_hg38.json"  // path to ExpansionHunter catalog, use 'default' for the default catalog

process REViewer {
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
    samtools sort -o ${sam}_sorted_realigned.bam 
    samtools index ${sam}_sorted_realigned.bam
    REViewer \
    --reads ${bamlet} \
    --vcf ${vcf} \
    --reference ${params.illumina_ref_fasta} \
    --catalog ${params.catalog} \
    --output-prefix ${sam}
    #--locus <Locus to analyze> \
    """
}


