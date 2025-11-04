params.caller = 'LongTR'
params.longtr_loci = "${projectDir}/catalogues/STRchive-disease-loci.hg38.longTR.bed"  

workflow run_longtr {
    take:
        sam_bam_ch
    main:
        longtr_results = sam_bam_ch |
            longtr
    emit:
        longtr_results  
}

process longtr {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "output/longtr/", mode: "copy", saveAs: { filename ->  filename.replaceAll("${sam}\\.", "${sam}_${type}.")}


    tag "${sam}_${type}"

    container 'quay.io/biocontainers/longtr:1.2--h077b44d_1'

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}.vcf.gz")

    script:
    """
    LongTR --bams ${bam} --fasta ${params.ref_fasta} --regions ${params.longtr_loci} --tr-vcf ${sam}.vcf.gz --min-mean-qual 10
    # --phased-bam
    """
}

