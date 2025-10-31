params.caller = 'ExpnansionHunter'
params.catalog = '../catalogues/scattr_hg38.tsv'

include { call } from './call.nf'

workflow run_scattr {
    take:
        sam_bam_ch
    main:
        scattr_results = sam_bam_ch |
            call
    emit:
        scattr_results  
}