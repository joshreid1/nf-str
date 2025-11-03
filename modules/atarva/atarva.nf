params.caller = 'atarva'
params.atarva_loci = "${projectDir}/catalogues/STRchive-disease-loci.hg38.atarva.bed"  

workflow run_atarva {
    take:
        sam_bam_ch
    main:
        trgt_results = sam_bam_ch |
            trgt
    emit:
        trgt_results  
}

process atarva {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "progress/atarva/", mode: "symlink", saveAs: { filename ->  filename.replaceAll("${sam}\\.", "${sam}_${type}.")}


    tag "${sam}_${type}"

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}.vcf")

    script:
    """
    atarva --fasta ${params.ref_fasta} \
           --bam ${bam} \
           --regions ${params.atarva_loci} \
           --vcf ${sam}.vcf
    """
}

