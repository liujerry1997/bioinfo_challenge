// tempus_bioinfo_challenge

// Set up other params
date = new Date().format("yyyy_MM_dd__hh_mm_ss")
params.publish_dir = "${params.publish_base}/${date}"

include { ANNOTATE_VCF } from './modules/annotate_vcf'

workflow {
    ANNOTATE_VCF(params.vcf_file_path, params.output_vcf_name)
}