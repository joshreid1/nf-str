params.caller = 'TRGT'
params.trgt_loci = "${projectDir}/catalogues/STRchive-disease-loci.hg38.TRGT.bed"  

workflow run_trgt {
    take:
        sam_bam_ch
    main:
        trgt_results = sam_bam_ch |
            trgt
    emit:
        trgt_results  
}

process trgt {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "output/trgt/", mode: "copy", saveAs: { filename ->  filename.replaceAll("${sam}\\.", "${sam}_${type}.")}


    tag "${sam}_${type}"

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}.vcf.gz"), path("${sam}.spanning.bam")

    script:
    """
    trgt genotype --genome ${params.ref_fasta} \
                  --repeats ${params.trgt_loci} \
                  --reads ${bam} \
                  --output-prefix ${sam}
    """
}

