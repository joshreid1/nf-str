params.caller = 'atarva'
params.atarva_loci = "${projectDir}/catalogues/STRchive-disease-loci.hg38.atarva.sorted.bed.gz"  

workflow run_atarva {
    take:
        sam_bam_ch
    main:
        atarva_results = sam_bam_ch |
            atarva
    emit:
        atarva_results  
}

process atarva {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "output/atarva/", mode: "copy", saveAs: { filename ->  filename.replaceAll("${sam}\\.", "${sam}_${type}.")}


    tag "${sam}_${type}"

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}.vcf")

    script:
    """
    micromamba run -n atarva atarva --fasta ${params.ref_fasta} \
           --bam ${bam} \
           --regions ${params.atarva_loci} \
           --vcf ${sam}.vcf
    """
}

