params.caller = 'Straglr'
params.straglr_loci = "${projectDir}/catalogues/STRchive-disease-loci.hg38.straglr.bed"  

workflow run_straglr {
    take:
        sam_bam_ch
    main:
        straglr_results = sam_bam_ch |
            straglr
    emit:
        straglr_results  
}

process straglr {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "progress/straglr/", mode: "symlink"
    tag { sam }

    //container '/vast/scratch/users/reid.j/nf-str-run/straglr_1.5.5--pyhdfd78af_0.sif'
    container 'quay.io/biocontainers/straglr:1.5.5--pyhdfd78af_0'

    input:
        tuple val(sam), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}.tsv"), path("${sam}.bed"), path("${sam}.vcf")

    script:
    """
        python /usr/local/bin/straglr.py ${bam} ${params.ref_fasta} ${sam} --loci ${params.straglr_loci} 
       
        #[--loci loci.bed] [--exclude skip_regions.bed] [--chroms chr] [--regions regions.bed] \
        #[--min_support N] [--min_ins_size N] [--min_str_len N] [--max_str_len N] [--nprocs N] \
        #[--genotype_in_size] [--max_num_clusters N] [--min_cluster_size N] [--working_dir] [--tmpdir] [--debug]
    """
}