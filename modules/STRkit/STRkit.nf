params.caller = 'STRkit'
params.strkit_loci = "${projectDir}/catalogues/STRkit_pathogenic_assoc.hg38.tsv"  

workflow run_strkit {
    take:
        sam_bam_ch
    main:
        strkit_results = sam_bam_ch |
            strkit
    emit:
        strkit_results  
}

process strkit {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "output/strkit/", mode: "copy"

    tag "${sam}_${type}"

    container 'ghcr.io/davidlougheed/strkit:latest'

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}_${type}.vcf")

    script:
    """
    strkit call \
    ${bam} --ref ${params.ref_fasta} --loci ${params.strkit_loci} \
    --vcf ${sam}_${type}.vcf --seed 123 \
    --hq --realign --no-tsv
    #--processes X --use-hp --incorporate-snvs path/to/dbsnp/00-common_all.vcf.gz 
    """
}

