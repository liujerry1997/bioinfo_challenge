// tempus_bioinfo_challenge

// Set up other params
date = new Date().format("yyyy_MM_dd__hh_mm_ss")
params.publish_dir = "${params.publish_base}/${date}"


process ANNOTATE_VCF {
    // Run python script to annotate variants

    publishDir "${params.publish_dir}", mode: "copy"
    cpus 1
    memory "2 GB"
    debug true
    
    input:
        path vcf_file_path
        val output_tsv_name
    output:
        path output_tsv_name, emit: txt
    script:
    """
    annotate_vcf.py --vcf_file_path ${vcf_file_path} --output_tsv_name ${output_tsv_name}
    """
}

workflow {
    ANNOTATE_VCF(params.vcf_file_path, params.output_tsv_name)
}