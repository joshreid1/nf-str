params.catalog = 'default'  // path to ExpansionHunter catalog, use 'default' for the default catalog

process call {
    cpus 8
	memory {'8 GB'}
	time '10 h'
    publishDir "output/scattr/", mode: "copy"
    tag { sam }

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), 
        path("${sam}.stats.json"), 
        path("${sam}.insert_distr.png"), 
        path("${sam}.depth_distr.png"),
        path("${sam}.bag.bam"),
        path("${sam}.defs.json"),
        path("${sam}.genotypes.json")

    script:
    """
    # Step 1: stats
    scattr ${sam} stats -@ 8 ${bam}

    # Step 2: extract
    scattr  extract -@ 8 ${bam} 

    # Step 3: define
    scattr ${sam} define ${params.catalog} ${params.ref_fasta}

    # Step 4: genotype
    scattr ${sam} genotype
    """
}