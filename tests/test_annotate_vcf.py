import filecmp
from bin.annotate_vcf import (
    annotate_variants,
    get_variant_ref_read_percentage,
    get_variant_allele_percentage,
    parse_helper,
    get_minor_allele_frequency,
)


def test_get_variant_ref_read_percentage():
    """
    A function to test the calculation of the percentage of reads supporting
    the variant versus those supporting reference reads.
    """
    assert get_variant_ref_read_percentage(total_coverage=10, reads_supporting_variant=5) == 100
    assert get_variant_ref_read_percentage(total_coverage=1, reads_supporting_variant=0) == 0
    assert get_variant_ref_read_percentage(total_coverage=1, reads_supporting_variant=1) == 100
    # No total variant coverage
    assert get_variant_ref_read_percentage(total_coverage=0, reads_supporting_variant=1) == None
    # Unexpected variant coverage
    assert get_variant_ref_read_percentage(total_coverage=10, reads_supporting_variant=-1) == None


def test_get_variant_allele_percentage():
    """
    A function to test the calculation of the percentage of reads supporting the variant.
    """
    assert get_variant_allele_percentage(10, 5) == 50
    assert get_variant_allele_percentage(1, 0) == 0
    assert get_variant_allele_percentage(1, 1) == 100
    assert get_variant_allele_percentage(0, 1) == None  # No total variant coverage
    assert get_variant_allele_percentage(10, -1) == None  # Unexpected variant coverage


def test_get_minor_allele_frequency():
    """
    A function to test the calculation of the minor allele frequency.
    """
    ref = "A"
    alt_list = ["T"]
    genotype_list = [[0, 0, False]]
    assert get_minor_allele_frequency(ref, alt_list, genotype_list) == 0

    ref = "A"
    alt_list = ["T", "C"]
    genotype_list = [[1, 0, False]]
    assert get_minor_allele_frequency(ref, alt_list, genotype_list) == 0

    ref = "A"
    alt_list = ["T", "C"]
    genotype_list = [[1, 0, False], [1, 2, False]]
    assert get_minor_allele_frequency(ref, alt_list, genotype_list) == 0.25

    ref = "A"
    alt_list = ["T", "C"]
    genotype_list = [[1, 0, False], [2, 2, False]]
    assert get_minor_allele_frequency(ref, alt_list, genotype_list) == 0.25

    # Missing genotype_list
    ref = "A"
    alt_list = ["T", "C"]
    genotype_list = None
    assert get_minor_allele_frequency(ref, alt_list, genotype_list) == None

    # Unexpected allele number
    ref = "A"
    alt_list = ["T", "C"]
    genotype_list = [[0, 3, False]]
    assert get_minor_allele_frequency(ref, alt_list, genotype_list) == None


def test_parse_helper():
    """
    A function to test the parsing of the variant information.
    """
    assert parse_helper(1) == 1
    assert parse_helper("one") == None  # Not an int
    assert parse_helper(1) == 1
    assert parse_helper((1, 1)) == 1
    assert parse_helper((1, 2)) == 1  # Take the integer digit of the average
    assert parse_helper((1, None)) == None  # Second value is not an int


def test_annotate_variants():
    """
    A function to test if annotate_variants function properly annotates variants.
    """
    vcf_file_path = "tests/data/short_vcf_data.txt"
    output_tsv_name = "tests/data/short_vcf_data_annotated.tsv"
    output_tsv_name = annotate_variants(vcf_file_path, output_tsv_name)

    expected_tsv = "tests/data/expected_output.tsv"
    assert filecmp.cmp(output_tsv_name, expected_tsv)
