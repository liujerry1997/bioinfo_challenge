
process ANNOTATE_VCF {
    // Run python script to annotate variants

    publishDir "${params.publish_dir}", mode: "copy"
    cpus 1
    memory "2 GB"
    debug true
    
    input:
        path vcf_file_path
        val output_vcf_name
    output:
        path output_vcf_name, emit: txt
    script:
    """
    annotate_vcf.py --vcf_file_path ${vcf_file_path} --output_vcf_name ${output_vcf_name}
    """
}

