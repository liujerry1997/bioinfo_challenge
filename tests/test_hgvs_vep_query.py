import unittest
from unittest.mock import patch
from bin.hgvs_vep_query import get_variant_vep_info


class TestGetVariantVepInfo(unittest.TestCase):
    @patch("bin.hgvs_vep_query.requests.get")
    def test_get_variant_vep_info(self, mock_get):
        mock_response = [
            {
                "allele_string": "G/C",
                "seq_region_name": "9",
                "assembly_name": "GRCh38",
                "start": 22125504,
                "strand": 1,
                "end": 22125504,
                "id": "9:g.22125504G>C",
                "input": "9:g.22125504G>C",
                "most_severe_consequence": "intron_variant",
                "transcript_consequences": [
                    {
                        "strand": 1,
                        "hgnc_id": "HGNC:34341",
                        "consequence_terms": ["downstream_gene_variant"],
                        "transcript_id": "ENST00000421632",
                        "impact": "MODIFIER",
                        "distance": 4815,
                        "variant_allele": "C",
                        "gene_symbol_source": "HGNC",
                        "biotype": "lncRNA",
                        "gene_symbol": "CDKN2B-AS1",
                        "gene_id": "ENSG00000240498",
                    },
                ],
            }
        ]

        mock_get.return_value.ok = True
        mock_get.return_value.json.return_value = mock_response

        hgvs_notation = "9:g.22125504G>C"
        vep_info = get_variant_vep_info(hgvs_notation)
        expected_vep_info = {
            "gene_id": "ENSG00000240498",
            "variation_type": "MODIFIER",
            "effect": "intron_variant",
        }
        self.assertEqual(vep_info, expected_vep_info)


if __name__ == "__main__":
    unittest.main()
