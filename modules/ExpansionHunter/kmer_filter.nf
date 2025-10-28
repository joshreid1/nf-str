params.catalog = "${projectDir}/catalogues/ExpansionHunter_hg38.json"  

process kmer_filter {
    
    publishDir "progress/eh5/kmer_filter/", mode: "symlink"
    tag { sam }

    input:
    tuple val(sam), path(bamlet_srt), path(bamlet_bai), path(vcf)

    output:
    tuple val(sam), path("${sam}_validated.vcf"), path("${sam}_validated.vcf.tbi")
    
    script:
    out_vcf = "${sam}_validated.vcf"
    out_vcf_tbi = "${sam}_validated.vcf.tbi"
    """
    micromamba run -n eh5_kmerfilter python kmer_filter.py --bam ${bamlet_srt} --vcf ${vcf} --catalog ${params.catalog} --auto --keep_lowdepth \
-o ./${sam}
    micromamba run -n base tabix -p vcf ${out_vcf}
    """
}