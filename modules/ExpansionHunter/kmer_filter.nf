params.catalog = "${projectDir}/catalogues/ExpansionHunter_hg38.json"  
process kmer_filter {
    
    publishDir "progress/eh5/kmer_filter/", mode: "symlink"
    
    input:
    tuple val(sam), path(bamlet_srt), path(bamlet_bai), path(vcf), path(catalog)
    
    output:
    tuple val(sam), path(out_vcf), path(out_vcf_tbi)
    
    script:
    out_vcf = "${sam}_validated.vcf"
    out_vcf_tbi = "${sam}_validated.vcf.tbi"
    """
    micromamba run -n eh5_kmerfilter python kmer_filter.py --bam ${bamlet_srt} --vcf ${vcf} --catalog ${catalog} --auto --keep_lowdepth \
-o ./${sam}
    micromamba run -n base tabix -p vcf ${out_vcf}
    """
}