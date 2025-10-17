params.catalog = 'default'  // path to ExpansionHunter catalog, use 'default' for the default catalog

process call {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "progress/eh5/", mode: "symlink"
    tag { sam }

    input:
        tuple val(sam), path(bam), path(bai), path(ref)

    output:
        tuple val(sam), path("${sam}_relaigned.bam"), path("${sam}.json"), path("${sam}.vcf")

    script:
    """
        ExpansionHunter \
		--reads ${bam} \
		--reference ${ref} \
		--variant-catalog ${params.catalog} \
		--output-prefix ${sam}
        """
}