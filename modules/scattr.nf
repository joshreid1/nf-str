params.caller = 'ExpnansionHunter'
params.catalog = '../catalogues/scattr_hg38.tsv'

include { call } from './scattr/call.nf'

workflow run_scattr {
    take:
        ref
        sam_bam_ch
    main:
        scattr_results = sam_bam_ch |
            combine(ref) |
            call
    emit:
        scattr_results  
}