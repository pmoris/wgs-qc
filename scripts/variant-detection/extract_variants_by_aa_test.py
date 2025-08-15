import unittest

target = __import__("extract_variants_by_aa")

# from .extract_variants import snpeff_annotation_2_dict

detect_markers_in_bcf_query = target.detect_markers_in_bcf_query
collect_annotated_samples = target.collect_annotated_samples
snpeff_annotation_2_dict = target.snpeff_annotation_2_dict

# Test if snpEff annotation is parsed correctly
marker_example = "G437A"
annotation_list_example = [
    "C|missense_variant|MODERATE|PF3D7_9999999|PF3D7_9999999|transcript|PF3D7_9999999.1|protein_coding|2/3|c.1310G>C|p.Gly437Ala|1614/2866|1310/2121|437/706||",
    "T|missense_variant|MODERATE|PF3D7_9999999|PF3D7_9999999|transcript|PF3D7_9999999.1|protein_coding|2/3|c.1310G>C|p.Gly437Ala|1614/2866|1310/2121|437/706||",
    "C|upstream_gene_variant|MODIFIER|PF3D7_0810700|PF3D7_0810700|transcript|PF3D7_0810700.1|protein_coding||c.-2320C>G|||||2235|",
    "C|upstream_gene_variant|MODIFIER|PF3D7_0811000|PF3D7_0811000|transcript|PF3D7_0811000.1|protein_coding||c.-4752G>C|||||3431|",
    "C|downstream_gene_variant|MODIFIER|PF3D7_0810600|PF3D7_0810600|transcript|PF3D7_0810600.1|protein_coding||c.*4889G>C|||||3711|",
    "C|downstream_gene_variant|MODIFIER|PF3D7_0810900|PF3D7_0810900|transcript|PF3D7_0810900.1|protein_coding||c.*1514C>G|||||1217|",
    "C|downstream_gene_variant|MODIFIER|PF3D7_0810900|PF3D7_0810900|transcript|PF3D7_0810900.2|protein_coding||c.*2182C>G|||||2182|",
]
mutation_allele_dict = snpeff_annotation_2_dict(annotation_list_example, marker_example)
print(mutation_allele_dict)
assert mutation_allele_dict == {"G437A": {"C": None, "T": None}}

# test if samples are detected correctly in bcftools query output
marker_of_interest_example = {
    "gene_id": "PF3D7_9999999",
    "gene_name": "some-gene",
    "marker": "G437A",
    "drug": "nonsense-drug",
}
query_example = """
Pf3D7_01_v3\t549685\tG\tC,T\tC|missense_variant|MODERATE|PF3D7_9999999|PF3D7_9999999|transcript|PF3D7_9999999.1|protein_coding|2/3|c.1310G>C|p.Gly437Ala|1614/2866|1310/2121|437/706||,T|missense_variant|MODERATE|PF3D7_9999999|PF3D7_9999999|transcript|PF3D7_9999999.1|protein_coding|2/3|c.1310G>C|p.Gly437Ala|1614/2866|1310/2121|437/706||,C|upstream_gene_variant|MODIFIER|PF3D7_0810700|PF3D7_0810700|transcript|PF3D7_0810700.1|protein_coding||c.-2320C>G|||||2235|,C|upstream_gene_variant|MODIFIER|PF3D7_0811000|PF3D7_0811000|transcript|PF3D7_0811000.1|protein_coding||c.-4752G>C|||||3431|,C|downstream_gene_variant|MODIFIER|PF3D7_0810600|PF3D7_0810600|transcript|PF3D7_0810600.1|protein_coding||c.*4889G>C|||||3711|,C|downstream_gene_variant|MODIFIER|PF3D7_0810900|PF3D7_0810900|transcript|PF3D7_0810900.1|protein_coding||c.*1514C>G|||||1217|,C|downstream_gene_variant|MODIFIER|PF3D7_0810900|PF3D7_0810900|transcript|PF3D7_0810900.2|protein_coding||c.*2182C>G|||||2182|\t106264-002-077=0/0;106264-002-078=0/1;106264-002-079=1/1;106264-002-080=0/2;106264-002-081=1/2;106264-002-082=2/2;
"""
detected_markers_per_sample = {}

detected_markers = detect_markers_in_bcf_query(
    query_example, marker_of_interest_example, detected_markers_per_sample
)

results_example = [
    {
        "sample": "106264-002-078",
        "mutation": "G437A",
        "marker_ALT": "C",
        "genotype": "0/1",
        "alt_number": 1,
        "gene_name": "some-gene",
        "gene_id": "PF3D7_9999999",
    },
    {
        "sample": "106264-002-079",
        "mutation": "G437A",
        "marker_ALT": "C",
        "genotype": "1/1",
        "alt_number": 1,
        "gene_name": "some-gene",
        "gene_id": "PF3D7_9999999",
    },
    {
        "sample": "106264-002-080",
        "mutation": "G437A",
        "marker_ALT": "T",
        "genotype": "0/2",
        "alt_number": 2,
        "gene_name": "some-gene",
        "gene_id": "PF3D7_9999999",
    },
    {
        "sample": "106264-002-081",
        "mutation": "G437A",
        "marker_ALT": "C",
        "genotype": "1/2",
        "alt_number": 1,
        "gene_name": "some-gene",
        "gene_id": "PF3D7_9999999",
    },
    {
        "sample": "106264-002-081",
        "mutation": "G437A",
        "marker_ALT": "T",
        "genotype": "1/2",
        "alt_number": 2,
        "gene_name": "some-gene",
        "gene_id": "PF3D7_9999999",
    },
    {
        "sample": "106264-002-082",
        "mutation": "G437A",
        "marker_ALT": "T",
        "genotype": "2/2",
        "alt_number": 2,
        "gene_name": "some-gene",
        "gene_id": "PF3D7_9999999",
    },
]


for sample in detected_markers:
    assert sample in results_example
assert len(results_example) == len(detected_markers)
