#!/usr/bin/env python3

import requests
import logging


def get_variant_vep_info(hgvs_notation: str):
    """
    Retrieves variant information from VEP API using the provided HGVS notation.

    Args::
        hgvs_notation (str): The HGVS notation of the variant to query.

    Returns:
        dict: A dictionary containing gene_id, variation_type, and effect of the variant.
        If the request fails, returns a dictionary with default unknown values.
    """
    url = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs_notation}"
    headers = {"Content-Type": "application/json"}
    try:
        response = requests.get(url, headers=headers)
        if response.ok:
            return parse_variant_vep_info(response.json())
        else:
            return {"gene_id": "Unknown", "variation_type": "Unknown", "effect": "Unknown"}

    except requests.RequestException as e:
        # Log any request exceptions that occur
        logging.error(f"Request failed for {hgvs_notation}: {e}")
        return {"gene_id": "Unknown", "variation_type": "Unknown", "effect": "Unknown"}


def parse_variant_vep_info(vep_json: list) -> dict:
    """
    A function that parses variant information from VEP JSON response.

    Args:
        vep_json (list): A list containing variant information in JSON format.

    Returns:
        dict: A dictionary containing gene_id, variation_type, and effect of the variant.
    """
    json_dict = vep_json[0]  # Use the first entry in the response
    consequence_dict = {
        key: value for key, value in json_dict.items() if key.endswith("consequences")
    }
    one_consequence = list(consequence_dict.keys())[0]  # Take the first key
    one_consequence_dict = consequence_dict.get(one_consequence)[
        0
    ]  # Take the first dict in the list
    gene_id = one_consequence_dict.get("gene_id", "Unknown")
    variation_type = one_consequence_dict.get("impact", "Unknown")
    effect = json_dict.get("most_severe_consequence", "Unknown")

    return {"gene_id": gene_id, "variation_type": variation_type, "effect": effect}
