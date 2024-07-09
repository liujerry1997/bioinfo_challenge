#!/usr/bin/env python3

import argparse
import csv
import logging
from statistics import mean
import cyvcf2
from cyvcf2 import VCF
from hgvs_vep_query import get_variant_vep_info


added_head_info_list = [
    {
        "ID": "PctSuppReads",
        "Description": "Percentage of reads supporting the variant",
        "Type": "Float",
        "Number": "1",
    },
    {
        "ID": "Gene",
        "Description": "The gene of the variant",
        "Type": "String",
        "Number": ".",
    },
    {
        "ID": "VariationType",
        "Description": "The type of variation",
        "Type": "String",
        "Number": ".",
    },
    {
        "ID": "VariationEffect",
        "Description": "The consequence/effect of variation",
        "Type": "String",
        "Number": ".",
    },
    {
        "ID": "MAF",
        "Description": "minor allele frequency",
        "Type": "Float",
        "Number": "1",
    },
]


def annotate_variants(vcf_file_path: str, output_tsv_name: str) -> str:
    """
    Annotates variants in the VCF file with additional information.
    Added relevant annotations to output tsv file.

    Args:
        vcf_file_path (str): The path to the input VCF file.
        output_tsv_name (str): The name of the output VCF file with annotations.

    Returns:
        output_tsv_name (str): The name of the output VCF file with annotations.
    """
    vcf_reader = VCF(vcf_file_path)
    for info_header in added_head_info_list:
        vcf_reader.add_info_to_header(info_header)

    with open(output_tsv_name, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, delimiter="\t")
        # Write the header row
        writer.writerow(
            [
                "CHROM",
                "POS",
                "REF",
                "ALT",
                "HGVS Notation",
                "total_coverage",
                "reads_supporting_variant",
                "Read Pct Varaint verses Reference",
                "Variant Allele Freq",
                "Gene ids",
                "Variation Types",
                "Effects",
                "Genotypes",
                "Minor Allele Freq",
            ]
        )
        for variant in vcf_reader:
            writer.writerow(process_one_variant(variant))

    return output_tsv_name


def process_one_variant(variant: cyvcf2.cyvcf2.Variant) -> list:
    """
    Process a single variant by calculating various metrics and querying VEP for additional information.

    Args:
        variant (cyvcf2.cyvcf2.Variant): The variant object to be processed.

    Returns:
        list: A list containing information about the processed variant including chromosomal position,
              reference and alternate alleles, HGVS notation, percentage of reads supporting the variant verses reference,
              variant allele frequency, gene IDs, variation types, effects, genotypes, and minor allele frequency.
    """
    total_coverage = parse_helper(variant.INFO.get("TC"))
    reads_supporting_variant = parse_helper(variant.INFO.get("TR"))
    chrom, pos, ref, alt_list = variant.CHROM, variant.POS, variant.REF, variant.ALT
    genotype_list = variant.genotypes

    # Calculate the percentage of reads supporting the variant and add to INFO field
    percentage_supporting_reads = get_variant_ref_read_percentage(
        total_coverage, reads_supporting_variant
    )
    variant_allele_freq = get_variant_allele_percentage(total_coverage, reads_supporting_variant)

    # Query the VEP for the variant, and add the information to the INFO field
    gene_ids, variation_types, effects = [], [], []
    for alt in alt_list:  # Extacted ALT is a list of strings.
        hgvs_notation = f"{chrom}:g.{pos}{ref}>{alt}"
        vep_info_dict = get_variant_vep_info(hgvs_notation)
        if vep_info_dict:
            gene_ids.append(vep_info_dict.get("gene_id"))
            variation_types.append(vep_info_dict.get("variation_type"))
            effects.append(vep_info_dict.get("effect"))
    minor_allele_freq = get_minor_allele_frequency(ref, alt_list, genotype_list)

    return [
        chrom,
        pos,
        ref,
        alt_list,
        hgvs_notation,
        total_coverage,
        reads_supporting_variant,
        percentage_supporting_reads,
        variant_allele_freq,
        ",".join(gene_ids),
        ",".join(variation_types),
        ",".join(effects),
        genotype_list,
        minor_allele_freq,
    ]


def get_variant_ref_read_percentage(total_coverage: int, reads_supporting_variant: int) -> float:
    """
    Calculate the percentage of reads supporting the variant versus those supporting reference reads.

    Args:
        total_coverage (int): The total coverage at the variant locus.
        reads_supporting_variant (int): The number of reads supporting the variant.

    Returns:
        float: The percentage of reads supporting the variant versus those supporting reference reads.
    """
    if total_coverage > 0 and reads_supporting_variant >= 0:
        reads_supporting_ref = total_coverage - reads_supporting_variant
        if total_coverage == reads_supporting_variant:
            return 100
        percentage_var_ref = (reads_supporting_variant / reads_supporting_ref) * 100
    else:
        percentage_var_ref = None
    return percentage_var_ref


def get_variant_allele_percentage(total_coverage: int, reads_supporting_variant: int) -> float:
    """
    Calculate the percentage of reads supporting the variant.

    Args:
        total_coverage (int): The total coverage at the variant locus.
        reads_supporting_variant (int): The number of reads supporting the variant.

    Returns:
        float: The percentage of reads supporting the variant.
    """
    if total_coverage > 0 and reads_supporting_variant >= 0:
        percentage_supporting_reads = (reads_supporting_variant / total_coverage) * 100
    else:
        percentage_supporting_reads = None
    return percentage_supporting_reads


def parse_helper(field_value):
    if isinstance(field_value, tuple):
        # If any value in the tuple is not an int, return None
        if any(not isinstance(val, int) for val in field_value):
            return None
        # Take the average int of the values
        return int(mean(field_value))
    elif isinstance(field_value, int):
        return field_value
    else:
        return None


def get_minor_allele_frequency(ref, alt_list: list[str], genotype_list: list) -> float:
    """
    Calculate the minor allele frequency based on the reference allele, list of alternate alleles, and list of genotypes.
    The given test vcf only has one sample per variant. The function is designed to handle multiple samples as well.

    Args:
        ref (str): The reference allele
        alt_list (list[str]): List of alternate alleles
        genotype_list (list): List of genotypes
    Returns:
        float: The minor allele frequency as a float
    """
    if genotype_list is None:
        return None
    allele_counts = {ref: 0}
    for alt in alt_list:
        allele_counts[alt] = 0

    for genotype in genotype_list:
        for allele in genotype:
            if allele > len(alt_list) or allele < 0:
                logging.error(f"Error: Unexpected allele number {allele} in genotype {genotype}")
                return None
            if isinstance(allele, bool):  # Skip the phase boolean
                continue
            if allele == 0:
                allele_counts[ref] += 1
            if allele > 0:
                allele_counts[alt_list[allele - 1]] += 1

    total_alleles = sum(allele_counts.values())
    if total_alleles == 0:
        return None

    frequencies = {allele: count / total_alleles for allele, count in allele_counts.items()}
    return min(frequencies.values())


def main():
    parser = argparse.ArgumentParser(
        description="Annotate VCF file with the percentage of reads supporting the variant."
    )
    parser.add_argument("--vcf_file_path", help="Path to the input VCF file")
    parser.add_argument("--output_tsv_name", help="Name of the output VCF file with annotations")

    args = parser.parse_args()
    annotate_variants(args.vcf_file_path, args.output_tsv_name)


if __name__ == "__main__":
    main()
