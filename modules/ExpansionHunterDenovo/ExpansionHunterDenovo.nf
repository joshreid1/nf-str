include { profile } from './profile.nf'
include { merge } from './merge.nf'

workflow run_expansion_hunter_denovo {
    take:
        sam_bam_ch
    main:
        // Run profile process
        profile_ch = sam_bam_ch | profile
        
        // Create manifest using collectFile operator
        manifest_ch = profile_ch
            .collectFile(name: 'ehdn_manifest.tsv', newLine: true) { sam, json_path ->
                "${sam}\tcontrol\t${json_path}"
            }
        
        // Pass manifest to merge process
        ehdn_results = merge(manifest_ch)

    emit:
        ehdn_results  
}