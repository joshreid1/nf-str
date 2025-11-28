process profile {
    cpus 1
	memory {'10 GB'}
	time '5 h'
    publishDir "output/ehdn/", mode: "copy"
    tag { sam }

    container 'quay.io/biocontainers/expansionhunterdenovo:0.9.0--h6ac36c1_11'

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}_${type}*json")

    script:
    """
        ExpansionHunterDenovo profile \
		--reads ${bam} \
		--reference ${params.illumina_ref_fasta} \
		--output-prefix ${sam}_${type} \
        #--min-anchor-mapq 50 \
		#--max-irr-mapq 40
    """
}

