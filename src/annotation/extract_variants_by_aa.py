# !TODO: create index if it does not exist

import argparse
import re
import subprocess
import warnings
from pathlib import Path

import pandas as pd


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
    rf"^{[''.join(AA_dict.keys())]}{{1}}\d+{[''.join(AA_dict.keys())]}{{1}}$"
)


def query_vcf(vcf_file, chrom, start, end):
    """Run bcftools query to extract variants within a specific region 'chrom:start-end'."""

    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\t[%SAMPLE=%GT;]\n' {vcf_file} -r {chrom}:{start}-{end}"

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

                gff_annotations[gene_id] = (chrom, start, end)

    return gff_annotations


def snpeff_annotation_2_dict(annotation_list, marker):
    """Returns a dictionary of mutations (should only be a single one)
    each of which contains a dictionary of ANN_Alleles, i.e. the
    alt alleles in the snpEff annotations of a particular location
    that correspond to the same marker mutation. The values of each
    alt allele is None, since this is set to the position of the ALT allele
    in the 5th VCF column (for matching with the genotype of each sample) later
    in collect_annotated_samples().

    E.g. (fictional example, codon is incorrect)
    Pf3D7_01_v3     30174   .       T       C,A     502.89  PASS    AC=2,2;AF=0.500,0.500;AN=4;DP=57;ExcessHet=0.0000;FS=0.000;MLEAC=4,6;MLEAF=1.00,1.00;MQ=46.83;QD=32.67;SOR=1.270;ANN=
    A|missense_variant|MODERATE|PF3D7_0100100|PF3D7_0100100|transcript|PF3D7_0100100.1|protein_coding|1/2|c.665T>A|p.Ile222Lys|665/6492|665/6492|222/2163||,
    T|missense_variant|MODERATE|PF3D7_0100100|PF3D7_0100100|transcript|PF3D7_0100100.1|protein_coding|1/2|c.665T>A|p.Ile222Lys|665/6492|665/6492|222/2163||
    C|missense_variant|MODERATE|PF3D7_0100100|PF3D7_0100100|transcript|PF3D7_0100100.1|protein_coding|1/2|c.665T>C|p.Ile222Thr|665/6492|665/6492|222/2163||      GT:AD:DP:GQ:PGT:PID:PL:PS       1|1:0,3,0:3:9:1|1:30161_AGT_A:135,9,0,135,9,135:30161   ./.     ./.     ./.:1,0,0:1:0:.:.:0,0,0,0,0,0   ./.     ./.:2,0,0:2:0:.:.:0,0,0,0,0,0        ./.     ./.:15,0,0:15:0:.:.:0,0,0,0,0,0 ./.     ./.     ./.:1,0,0:1:0:.:.:0,0,0,0,0,0   2|2:0,0,8:8:24:1|1:30173_A_G:360,360,360,24,24,0:30173  ./.     ./.     ./.     ./.

    would results in the dictionary:
    {
    I222K : {
        A : None,
        T : None,
        }
    }


    Args:
        annotation_list (_type_): _description_
        marker (_type_): _description_

    Returns:
        _type_: _description_
    """
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

        # skip annotations that do not contain a protein-level variant annotation
        if not ANN_HGVS_p:
            continue

        # convert three letter code into 1 letter codes to match variants of interest list
        mutation = ANN_HGVS_p
        for A, AAA in AA_dict.items():
            mutation = mutation.replace(AAA, A)
        mutation = mutation.replace("p.", "")

        # ignore other types of mutations like delins or fs, since they do not occur in the variants list
        if not r.match(mutation):
            warnings.warn(f"{ANN_HGVS_p} does not match expected format...")
            continue
        # TODO: how to handle p.Leu261delinsTyrIle - disruptive_inframe_insertion ?
        # p.Met74fs does not match expected format - frameshift

        # skip if annotated mutation does not match marker in variants of interest list
        # TODO: this can probably replace the above two conditionals
        if not mutation == marker:
            continue

        # store all the different annotations that cause the same mutation under their ALT allele
        if mutation not in mutation_allele_dict:
            mutation_allele_dict[mutation] = {ANN_Allele: None}
        # check if ANN_Allele is present twice in snpEff annotations for a single position, should not be possible...
        elif ANN_Allele in mutation_allele_dict[mutation].keys():
            raise ValueError(
                f"ANN_Allele {ANN_Allele} for mutation {mutation} was found twice in snpEff annotation. Should not be possible! \n {mutation_allele_dict}"
            )
        mutation_allele_dict[mutation][ANN_Allele] = None

        # check if there are multiple alleles leading to the same marker mutation, since these are noteworthy
        if len(mutation_allele_dict[mutation].keys()) > 1:
            print(
                f"Found multiple annotations for the same mutation! {marker}"
            )  # {mutation_allele_dict[mutation]}

        # TODO: check if multiple annotations with the same alt ANN_Allele can occur. If so, needs to be taken into account so that they do not overwrite each other.
        # ! TODO: write test input file to check different conditions like this

    return mutation_allele_dict


def collect_annotated_samples(
    mutation_allele_dict, samples, alt, marker_of_interest, detected_markers_per_sample
):
    """Find all samples that contain any of the various annotations for a specific mutation
    and return them as a list of dictionaries (plus update global dictionary for
    simplified form for printing).
    For a sample with a genotype where both alleles can cause the same mutation, both will
    be added as individual dictionaries (= rows in eventual dataframe).

    Args:
        mutation_allele_dict (_type_): A dictionary containing mutations (as keys, should always
        only be a single one) and different ALT alleles that correspond to it.
        samples (_type_): _description_

    Returns:
        _type_: _description_
    """
    detected_markers_list = []

    # iterate over sample genotypes for current row of VCF ( = bcftools query output)
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

        alleles = genotype.split("/") if "/" in genotype else genotype.split("|")

        print(f"Found marker in {sample_name} with genotype {genotype}")
        #! print(line)

        # iterate over the two alleles (assumes diploid calling)
        for allele in set(alleles):  # set avoids repeat check for homozygous genotypes

            # skip uncalled alleles or reference alleles
            if allele == "." or allele == "0":
                continue

            # warn if allele is not an integer (= position of allele in ALT column)
            try:
                allele = int(allele)
            except ValueError:
                warnings.warn(
                    f"Allele {allele} in genotype {genotype} could not be assigned to marker. Expected an integer to represent the position of the allele in the VCF ALT field."
                )
                continue

            # iterate over the different annotated alleles that cause the current mutation
            # {'G437A': {'C': None, 'T': None}}
            # should only ever contain a single mutation with 1 or more ANN_Alleles
            for mutation, ann_allele_dict in mutation_allele_dict.items():

                # check which of the ALT alleles is present in the current sample
                for ann_allele, alt_number in ann_allele_dict.items():
                    # skip if alt does not match allele as labeled in genotype
                    # e.g., 0/1, where 1 = alt_number taken from the position of the ANN_Allele,
                    # matching the current annotation/mutation combination,
                    # in the ALT column
                    if alt_number != allele:
                        # TODO: remove warning: this will happen every time there is a heterozygote with multiple alleles causing the same mutation, because we are checking each allele and annotation one by one, e.g. for 1/2, 1 will never match the annotation for 2 and vice versa.
                        warnings.warn(
                            f"Allele {allele} in genotype {genotype} could not be found in annotation {ann_allele}-{alt_number} ( {mutation_allele_dict}. Expected {alt_number} )."
                        )
                        continue

                    detected_marker = {
                        "sample": sample_name,
                        "mutation": mutation,
                        "marker_ALT": ann_allele,
                        "genotype": genotype,
                        "alt_number": alt_number,
                        "gene_name": marker_of_interest["gene_name"],
                        "gene_id": marker_of_interest["gene_id"],
                    }
                    assert ann_allele in alt
                    assert alt.split(",").index(ann_allele) + 1 == alt_number

                    detected_markers_list.append(detected_marker)

                    # create simple per sample summary
                    if sample_name not in detected_markers_per_sample:
                        detected_markers_per_sample[sample_name] = []
                    detected_markers_per_sample[sample_name].append(
                        f"{mutation} ({marker_of_interest["gene_name"]})"
                    )

    #!                            print(detected_markers)

    return detected_markers_list


def detect_markers_in_bcf_query(
    bcf_query_result, marker_of_interest, detected_markers_per_sample
):
    """Search for given marker in bcftools query output and return a list
    of dictionaries, each containing a single sample that has one of the
    various genotypes/annotations that correspond to the marker mutation,
    (plus store them in the global dictionary for easy printing).

    Args:
        bcf_query_result (_type_): _description_
        marker (_type_): _description_
        samples (_type_): _description_

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
    """

    # loop through lines of bcf query
    for line in bcf_query_result.strip().split("\n"):
        fields = line.split("\t")

        # check if vcf output format matches expectations
        if len(fields) != 6:
            raise ValueError(
                f"bcftools query output is not in the expected format: {line}. Expected: 6 tab-separated values (with the last value being a semicolon-separated list of sample genotypes)."
            )

        chrom, pos, ref, alt, annotation, samples = fields

        # variant can have multiple snpeff annotations
        annotation_list = annotation.split(",")

        # search for marker in snpeff annotations
        # multiple annotations can have the same mutation, so store them in a dictionary
        # {'G437A': {'C': None, 'T': None}}
        # should only ever contain a single mutation with 1 or more ANN_Alleles
        mutation_allele_dict = snpeff_annotation_2_dict(
            annotation_list, marker_of_interest["marker"]
        )

        # skip to next row if marker in variants of interest list could not be found in
        # annotations of current row of bcftools query output
        if not mutation_allele_dict:
            continue

        # handle multi-allelic sites by storing the position of the ANN alt in the main VCF ALT column
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

        print(
            f"Found annotation for {marker_of_interest["marker"]} in VCF file - querried on gene_id {marker_of_interest["gene_id"]} - current position = {pos}...\nLooking for matching genotype in samples..."
        )

        if len(mutation_allele_dict.keys()) > 1:
            raise ValueError(
                "Found multiple ANN mutations in a single query, should not be possible."
            )

        # store samples that match the found marker annotation
        detected_markers = collect_annotated_samples(
            mutation_allele_dict,
            samples,
            alt,
            marker_of_interest,
            detected_markers_per_sample,
        )

        return detected_markers


def parse_markers_file(marker_file):
    markers = pd.read_csv(marker_file, sep=",", skip_blank_lines=True, comment="#")
    markers = markers.dropna(
        subset=["gene_id", "gene_name", "mutation", "drug"]
    ).reset_index(drop=True)
    return markers


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract variants from a VCF file based on a list of markers of interest."
    )
    parser.add_argument(
        "-v",
        "--vcf",
        type=str,
        help="Input VCF file in .gz format. Requires index.",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--markers",
        type=str,
        help="CSV file with gene_id,gene_name and mutation (three-letter HGVS.p syntax) columns",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--gff",
        help="GFF annotation file, used for retrieving genomic coordinates of genes containing markers.",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path where to store output CSV file.",
        required=False,
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing output CSV file."
    )
    # parser.add_argument(
    #     "--expand", type=int, default=0, help="Expand region by ±N bp (default: 0)"
    # )
    args = parser.parse_args()

    # check input and output files
    vcf_file = Path(args.vcf)
    marker_file = Path(args.markers)
    gff_file = Path(args.gff)
    if not vcf_file.exists():
        raise FileNotFoundError(f"Input VCF file could not found: {vcf_file}")
    if not marker_file.exists():
        raise FileNotFoundError(f"Input markers CSV could not found: {marker_file}")
    if not gff_file.exists():
        raise FileNotFoundError(f"Input GFF file could not found: {gff_file}")
    output_file = Path(args.output) if args.output else None
    if output_file and output_file.exists() and not args.overwrite:
        warnings.warn(
            f"Output file already exists. Use --overwrite to force re-generation."
        )
        exit(1)

    # parse gff
    gff_annotations = parse_gff(gff_file)

    # parse variants/markers/mutations
    markers = parse_markers_file(marker_file)

    # create list (and dictionary) to store detected marker/sample data
    detected_markers = []
    detected_markers_per_sample = {}

    # TODO: add to function
    # iterate over each marker in marker list and search for its presence in the VCF file
    for _, marker_row in markers.iterrows():

        # gene_id, gene_name, marker, drug = (
        #     marker_row["gene_id"],
        #     marker_row["gene_name"],
        #     marker_row["mutation"],
        #     marker_row["drug"],
        #     # marker_row["info"],
        #     # TODO: add additional info depending on what is present in variants of interest file
        # )
        marker_of_interest = {
            "gene_id": marker_row["gene_id"],
            "gene_name": marker_row["gene_name"],
            "marker": marker_row["mutation"],
            "drug": marker_row["drug"],
        }

        # find marker gene start and end coordinates
        chrom, start, end = gff_annotations[marker_of_interest["gene_id"]]

        # extract vcf entries
        print(
            f"Querying VCF file for marker mutation {marker_of_interest["marker"]}..."
        )
        result = query_vcf(vcf_file, chrom, start, end)

        # skip empty output lines
        if not result.stdout:
            warnings.warn(
                f"bcftools query output was empty! {result.stdout} \n Search query region: '-r {chrom}:{start}-{end}'"
            )
            continue

        # detect markers in bcf query output
        marker_results = detect_markers_in_bcf_query(
            result.stdout, marker_of_interest, detected_markers_per_sample
        )
        if marker_results:
            # add additional column info back
            for result in marker_results:
                result["drug"] = marker_of_interest["drug"]
            # store in output list for conversion to dataframe
            detected_markers.extend(marker_results)

    # Convert results to DataFrame and save to CSV
    df_output = pd.DataFrame(detected_markers)
    df_output.to_csv(output_file, index=False) if output_file else None

    # TODO: also output wide instead of long format?

    # Output results
    print("\nDetected Markers Per Sample:")
    for sample, markers in detected_markers_per_sample.items():
        print(f"{sample}: {', '.join(markers)}")
