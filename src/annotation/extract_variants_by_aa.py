import subprocess
import warnings
import re
from pathlib import Path

import pandas as pd


# Function to run bcftools and extract nearby variants
def query_vcf(vcf_file, chrom, start, end):
    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%SAMPLE=%GT;]\n' {vcf_file} -r {chrom}:{start}-{end}"

    # print(cmd)

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error running bcftools query: {result.stderr}")

    # return result.stdout.strip().split("\n") if result.stdout else []
    return result


def parse_gff(gff_file):
    """Parses a GFF file to extract gene regions."""
    gff_annotations = {}

    with gff_file.open("r") as f:
        for line in f:
            # skip comment lines
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            # filter on CDS features
            if fields[2] == "protein_coding_gene":

                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]

                # split attribute field
                attributes = {
                    k: v
                    for k, v in (x.split("=") for x in fields[8].split(";") if "=" in x)
                }

                # extract gene id
                gene_id = attributes.get("ID", "").split(":")[-1]

                # add to dictionary if not yet present
                if gene_id in gff_annotations:
                    warnings.warn(
                        f"Duplicate entries found for protein coding gene entry: {line}"
                    )
                    # import ipdb

                    # ipdb.set_trace()

                    # gene_cds[gene_id] = {"chrom": chrom, "cds": [], "strand": strand}
                gff_annotations[gene_id] = (chrom, start, end)

    return gff_annotations


AA_dict = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}

r = re.compile(
    # rf"^{[''.join(AA_dict.keys())]}{{1}}\d{{1,5}}{[''.join(AA_dict.keys())]}{{1}}$"
    rf"^{[''.join(AA_dict.keys())]}{{1}}\d+{[''.join(AA_dict.keys())]}{{1}}$"
)

# Input files
vcf_file = Path(
    "/home/pmoris/itg/projects/summit/all-batches/batch003-pf/combined.filtered.ann.vcf.gz"
)  # Replace with your actual VCF file
marker_file = Path(
    "/home/pmoris/itg/projects/summit/molecular-markers/pf-resistance-simple.csv"
)  # Your list of known markers
gff_file = Path(
    "/home/pmoris/itg/projects/summit/all-batches/batch003-pf/PlasmoDB-68_Pfalciparum3D7.gff"
)

# gff_file = Path(args.gff_file)
# marker_file = Path(args.marker_file)
# bed_file = Path(args.bed_file)

# parse gff
gff_annotations = parse_gff(gff_file)

# Load marker list from TSV file
markers = pd.read_csv(marker_file, sep=",", skip_blank_lines=True, comment="#")
markers = markers.dropna(
    subset=["gene_id", "gene_name", "mutation", "drug"]
).reset_index(drop=True)

detected_markers = {
    "sample": [],
    "mutation": [],
    "marker_ALT": [],
    "genotype": [],
    "gene_name": [],
    "gene_id": [],
    "alt": [],
}

detected_markers_per_sample = {}

# iterate over each marker
for index, row in markers.iterrows():
    gene_id, gene_name, marker = (
        row["gene_id"],
        row["gene_name"],
        row["mutation"],
        # row["drug"],
        # row["info"],
    )

    print(f"Querying VCF file for marker mutation {marker}...")

    # find marker gene start and end coordinates
    chrom, start, end = gff_annotations[gene_id]

    # extract vcf entries
    result = query_vcf(vcf_file, chrom, start, end)

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            fields = line.split("\t")
            if len(fields) != 6:
                raise ValueError(
                    f"bcftools query output is not in the expected format: {line}. Expected: 5 tab-separated values (with the last value being a semicolon-separated list of sample genotypes)."
                )

            chrom, pos, ref, alt, annotation, samples = fields

            # variant can have multiple snpeff annotations
            annotation_list = annotation.split(",")

            # different annotations can have the same mutation
            mutation_allele_dict = {}

            for ann in annotation_list:
                (
                    ANN_Allele,
                    ANN_Annotation,
                    ANN_Annotation_Impact,
                    ANN_Gene_Name,
                    ANN_Gene_ID,
                    ANN_Feature_Type,
                    ANN_Feature_ID,
                    ANN_Transcript_BioType,
                    ANN_Rank,
                    ANN_HGVS_c,
                    ANN_HGVS_p,
                    ANN_cDNA_pos_ANN_cDNA_length,
                    ANN_CDS_pos_ANN_CDS_length,
                    ANN_AA_pos_ANN_AA_length,
                    ANN_Distance,
                    ANN_error,
                ) = ann.split("|")

                # TODO: handle other types of annotated mutations

                if not ANN_HGVS_p:
                    continue

                mutation = ANN_HGVS_p
                for A, AAA in AA_dict.items():
                    mutation = mutation.replace(AAA, A)
                mutation = mutation.replace("p.", "")

                # TODO: how to handle p.Leu261delinsTyrIle - disruptive_inframe_insertion ?
                # p.Met74fs does not match expected format - frameshift

                if not r.match(mutation):
                    warnings.warn(f"{ANN_HGVS_p} does not match expected format...")

                # skip if annotated mutation does not match marker in variants of interest list
                if not mutation == marker:
                    continue

                # TODO: check which annotation cause the same mutation
                # do we need to keep track of different mutations too?
                # also add check for current mutation in marker list loop
                #! print("reached here", ann, mutation, marker)

                # store different annotations that cause the same mutation
                if mutation not in mutation_allele_dict:
                    mutation_allele_dict[mutation] = {ANN_Allele: None}
                mutation_allele_dict[mutation][ANN_Allele] = None

                # TODO: check if multiple annotations with the same alt ANN_Allele can occur. If so, needs to be taken into account so that they do not overwrite each other.

            # handle multi-allelic sites
            alt_list = alt.split(",")
            for mutation, ann_allele_dict in mutation_allele_dict.items():
                for ann_allele in ann_allele_dict.keys():
                    try:
                        alt_number = alt_list.index(ann_allele)
                        mutation_allele_dict[mutation][ann_allele] = (
                            alt_number + 1
                        )  # 0 = REF, 1 = first ALT, 2 = second ALT, etc.
                    except ValueError:
                        warnings.warn(
                            f"\nCould not find allele in VCF ALT entry, despite it occuring in the SnpEff ANN output:\n{line}\n"
                        )

                # alt_number = alt_list.index(ANN_Allele[0]) if ANN_Allele[0] in alt_list else None
                # if alt_number: # ! does not behave as expected, since if 0 defaults to false
                # mutation_allele_dict[mutation].append(alt_number)
                # else:
                #     warnings.warn(
                #         f"\nCould not find allele in VCF ALT entry, despite it occuring in the SnpEff ANN output:\n{line}\n"
                #     )

            # skip if marker in variants of interest list could not be found in annotations of current row of bcftools query output
            if not marker in mutation_allele_dict.keys():
                # warnings.warn(
                #     f"Could not find {marker} in VCF file (querried on gene_id {gene_id})..."
                # )
                # print(line)
                # print("dict", mutation_allele_dict)
                # print("marker", marker)
                # exit(1)
                continue

            print(
                f"Found annotation for {marker} in VCF file - querried on gene_id {gene_id} - current position = {pos}..."
            )

            if len(mutation_allele_dict.keys()) > 1:
                raise ValueError(
                    "Found multiple ANN mutations in a single query, should not be possible."
                )

            for sample in samples.split(";"):

                # skip empty end of line
                if not sample:
                    #! print(samples)
                    # warnings.warn(f"Found empty sample column in VCF entry: {line}")
                    continue

                sample_name, genotype = sample.split("=")

                # print(f"Processing sample {sample_name}")

                # consider only samples with a non-reference allele (0/0)
                if genotype == "0/0" or genotype == "./.":
                    # print(f"No marker found in {sample_name}")
                    continue

                alleles = (
                    genotype.split("/") if "/" in genotype else genotype.split("|")
                )

                print(f"Found marker in {sample_name} with genotype {genotype}")
                #! print(line)

                # iterate over the two alleles (assumes diploid calling)
                # do not repeat for homozygous genotypes
                for allele in set(alleles):

                    # # skip if previous allele was the same
                    # if i == 1 and alleles[0] == alleles[1]:
                    #     continue

                    # skip uncalled alleles or reference alleles
                    if allele == "." or allele == "0":
                        continue

                    # warn if allele is not an integer (= position of allele in ALT column)
                    try:
                        allele = int(allele)
                    except ValueError:
                        warnings.warn(
                            f"Allele {allele} in genotype {genotype} could not be assigned to marker."
                        )
                        continue

                    # iterate over the different annotated alleles that cause the current mutation
                    for mutation, ann_allele_dict in mutation_allele_dict.items():

                        # check which ALT allele is present in the sample
                        for ann_allele, alt_number in ann_allele_dict.items():

                            #!                          print("checking allele", allele)
                            #!                          print(mutation, ann_allele)
                            #!                          print(alleles)
                            #!                          print(genotype)
                            #!                          print(alt_number)
                            #!                          print(line)
                            #!                          print("marker", marker)
                            #!                          print("dict", mutation, ann_allele_dict)

                            # skip if alt does not match allele as labeled in genotype (e.g., 0/1, where 1 = alt_number taken from the position of the ANN_Allele, matching the current annotation/mutation combination, in the ALT column
                            if alt_number != allele:
                                warnings.warn(
                                    f"Allele {allele} in genotype {genotype} could not found in annotations."
                                )
                                exit(1)
                                continue

                            detected_markers["sample"].append(sample_name)
                            detected_markers["mutation"].append(mutation)
                            detected_markers["marker_ALT"].append(ann_allele)
                            detected_markers["genotype"].append(genotype)
                            detected_markers["gene_name"].append(gene_name)
                            detected_markers["gene_id"].append(gene_id)
                            detected_markers["alt"].append(
                                alt_number
                            )  # TODO: check if position is correct
                            assert ANN_Allele in alt  # TODO
                            assert alt.index(ann_allele) + 1 == alt_number

                            # create simple per sample summary
                            if sample_name not in detected_markers_per_sample:
                                detected_markers_per_sample[sample_name] = []
                            detected_markers_per_sample[sample_name].append(
                                f"{mutation} ({gene_name})"
                            )

#!                            print(detected_markers)

# for allele in mutation_allele_dict.values():

# samples_genotypes = {}
# for sample in samples.split(";"):
#     if sample:
#         sample_name, genotype = sample.split("=")
#         samples_genotypes[sample_name] = genotype

#         # if "1" in genotype:
#         #     if sample_name not in sample_markers:
#         #         sample_markers[sample_name] = []
#         #     sample_markers[sample_name].append(f"{marker} ({gene_name})")

# # Consider only samples with a non-reference allele (e.g., 0/1 or 1/1)
# for sample_name, genotype in samples_genotypes.items():
#     if "1" in genotype:
#         if sample_name not in sample_markers:
#             sample_markers[sample_name] = []
#         sample_markers[sample_name].append(f"{marker} ({gene_name})")


# Convert results to DataFrame and save to CSV
df_output = pd.DataFrame(detected_markers)

# Output results
print("\nDetected Markers Per Sample:")
for sample, markers in detected_markers_per_sample.items():
    print(f"{sample}: {', '.join(markers)}")


# #####

# # Input files
# vcf_file = "/home/pmoris/itg/projects/summit/all-batches/batch003-pf/combined.filtered.ann.vcf.gz"  # Replace with your actual VCF file
# marker_file = "/home/pmoris/itg/projects/summit/all-batches/variants_pf.csv"  # Your list of known markers

# # Load marker list from TSV file
# markers = pd.read_csv(marker_file, sep=",", skip_blank_lines=True)
# markers = markers.dropna(subset=["CHROM", "POS", "Marker", "Gene_Name"]).reset_index(
#     drop=True
# )
# markers["POS"] = markers["POS"].astype(int)

# # Dictionary to store found markers per sample
# sample_markers = {}

# # Iterate over each marker
# for index, row in markers.iterrows():
#     chrom, pos, gene, gene_name, marker = (
#         row["gene_id"],
#         row["gene_name"],
#         row["mutation"],
#         # row["drug"],
#         # row["info"],
#     )

#     # Use bcftools to check if this variant is present in the VCF
#     cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE=%GT;]\n' {vcf_file} -r {chrom}:{pos}"
#     result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
#     if result.returncode != 0:
#         raise RuntimeError(f"Error running bcftools query: {result.stderr}")

#     if result.stdout:
#         for line in result.stdout.strip().split("\n"):
#             fields = line.split("\t")
#             if len(fields) != 5:
#                 raise ValueError(
#                     f"bcftools query output is not in the expected format: {line}. Expected: 5 tab-separated values (with the last value being a semicolon-separated list of sample genotypes)."
#                 )

#             chrom, pos, ref, alt, samples = fields

#             samples_genotypes = {}
#             for sample in samples.split(";"):
#                 if sample:
#                     sample_name, genotype = sample.split("=")
#                     samples_genotypes[sample_name] = genotype

#                     # if "1" in genotype:
#                     #     if sample_name not in sample_markers:
#                     #         sample_markers[sample_name] = []
#                     #     sample_markers[sample_name].append(f"{marker} ({gene_name})")

#             # Consider only samples with a non-reference allele (e.g., 0/1 or 1/1)
#             for sample_name, genotype in samples_genotypes.items():
#                 if "1" in genotype:
#                     if sample_name not in sample_markers:
#                         sample_markers[sample_name] = []
#                     sample_markers[sample_name].append(f"{marker} ({gene_name})")

# # Output results
# print("\nDetected Markers Per Sample:")
# for sample, markers in sample_markers.items():
#     print(f"{sample}: {', '.join(markers)}")
