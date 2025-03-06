import subprocess

import pandas as pd

# Input files
vcf_file = "/home/pmoris/itg/projects/summit/all-batches/batch003-pf/combined.filtered.ann.vcf.gz"  # Replace with your actual VCF file
marker_file = "/home/pmoris/itg/projects/summit/all-batches/variants_pf.csv"  # Your list of known markers

# Load marker list from TSV file
markers = pd.read_csv(marker_file, sep=",")

# Dictionary to store found markers per sample
sample_markers = {}

# Iterate over each marker
for index, row in markers.iterrows():
    chrom, pos, gene, gene_name, marker = (
        row["CHROM"],
        row["POS"],
        row["Gene_ID"],
        row["Gene_Name"],
        row["Marker"],
    )

    # Use bcftools to check if this variant is present in the VCF
    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE=%GT;]\n' {vcf_file} -r {chrom}:{pos}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Error running bcftools query: {result.stderr}")

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            fields = line.split("\t")
            if len(fields) != 5:
                raise ValueError(
                    f"bcftools query output is not in the expected format: {line}. Expected: 5 tab-separated values (with the last value being a semicolon-separated list of sample genotypes)."
                )

            chrom, pos, ref, alt, samples = fields

            samples_genotypes = {}
            for sample in samples.split(";"):
                if sample:
                    sample_name, genotype = sample.split("=")
                    samples_genotypes[sample_name] = genotype

                    # if "1" in genotype:
                    #     if sample_name not in sample_markers:
                    #         sample_markers[sample_name] = []
                    #     sample_markers[sample_name].append(f"{marker} ({gene_name})")

            # Consider only samples with a non-reference allele (e.g., 0/1 or 1/1)
            for sample_name, genotype in samples_genotypes.items():
                if "1" in genotype:
                    if sample_name not in sample_markers:
                        sample_markers[sample_name] = []
                    sample_markers[sample_name].append(f"{marker} ({gene_name})")

# Output results
print("\nDetected Markers Per Sample:")
for sample, markers in sample_markers.items():
    print(f"{sample}: {', '.join(markers)}")
