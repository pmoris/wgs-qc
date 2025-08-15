# TODO: change code so that it searches for all relevant files in a directory instead of needing a list of sample names
# needs to be aware that samples with similar naming structure need to be grouped

import argparse
import csv
import json
import re
from pathlib import Path

import pandas as pd


def extract_flagstat_data(flagstat_dir, sample_names):
    sample_flagstats = {}
    flagstat_dir = Path(flagstat_dir)

    for sample in sample_names:
        flagstat_file = flagstat_dir / f"{sample}.sort.markdup.bam.flagstat"

        if flagstat_file.exists():
            with flagstat_file.open("r") as f:
                lines = f.readlines()

                total_reads = int(lines[0].split("+")[0].strip())

                mapped_reads, mapped_percent = map(
                    str.strip,
                    re.findall(r"(\d+) \+ \d+ mapped \((\d+\.\d+)%", lines[6])[0],
                )

                primary_mapped_reads, primary_mapped_percent = map(
                    str.strip,
                    re.findall(r"(\d+) \+ \d+ primary mapped \((\d+\.\d+)%", lines[7])[
                        0
                    ],
                )

                properly_paired_reads, properly_paired_percent = map(
                    str.strip,
                    re.findall(
                        r"(\d+) \+ \d+ properly paired \((\d+\.\d+)%", lines[11]
                    )[0],
                )

                sample_flagstats[sample] = {
                    "Flagstat Total Reads": total_reads,
                    "Flagstat Mapped Reads": int(mapped_reads),
                    "Flagstat Mapped Percent": float(mapped_percent),
                    "Flagstat Primary Mapped Reads": int(primary_mapped_reads),
                    "Flagstat Primary Mapped Percent": float(primary_mapped_percent),
                    "Flagstat Properly Paired Reads": int(properly_paired_reads),
                    "Flagstat Properly Paired Percent": float(properly_paired_percent),
                }

    # for sample in sample_names:
    # flagstat_file = flagstat_dir / f"{sample}.sort.markdup.bam.flagstat"
    #
    # if flagstat_file.exists():
    # with flagstat_file.open('r') as f:
    # lines = f.readlines()
    #
    # total_reads = int(lines[0].split('+')[0].strip())
    # mapped_reads, mapped_percent = 0, 0.0
    # primary_mapped_reads, primary_mapped_percent = 0, 0.0
    # properly_paired_reads, properly_paired_percent = 0, 0.0
    #
    # for line in lines:
    # if " mapped (" in line:
    # mapped_reads, mapped_percent = map(str.strip, re.findall(r'(\d+) \+ 0 mapped \((\d+\.\d+)%', line)[0])
    # elif " primary mapped (" in line:
    # primary_mapped_reads, primary_mapped_percent = map(str.strip, re.findall(r'(\d+) \+ 0 primary mapped \((\d+\.\d+)%', line)[0])
    # elif " properly paired (" in line:
    # properly_paired_reads, properly_paired_percent = map(str.strip, re.findall(r'(\d+) \+ 0 properly paired \((\d+\.\d+)%', line)[0])
    #
    # sample_flagstats[sample] = {
    # 'total_reads': total_reads,
    # 'mapped_reads': int(mapped_reads),
    # 'mapped_percent': float(mapped_percent),
    # 'primary_mapped_reads': int(primary_mapped_reads),
    # 'primary_mapped_percent': float(primary_mapped_percent),
    # 'properly_paired_reads': int(properly_paired_reads),
    # 'properly_paired_percent': float(properly_paired_percent),
    # }

    return sample_flagstats


def extract_total_reads(json_dir, sample_names):
    sample_totals = {}
    json_dir = Path(json_dir)

    for sample in sample_names:
        sample_totals[sample] = {"before_filtering": 0, "after_filtering": 0}

        # Find all matching files for the sample
        json_files = list(json_dir.glob(f"*{sample}_*_L00*.trim.json"))

        for json_file in json_files:
            with json_file.open("r") as f:
                data = json.load(f)
                before_reads = data["summary"]["before_filtering"]["total_reads"]
                after_reads = data["summary"]["after_filtering"]["total_reads"]

                sample_totals[sample]["before_filtering"] += before_reads
                sample_totals[sample]["after_filtering"] += after_reads

    return sample_totals


def extract_screen_data(screen_dir, sample_names):
    sample_screens = {}
    screen_dir = Path(screen_dir)

    for sample in sample_names:
        sample_screens[sample] = {}

        screen_files = list(screen_dir.glob(f"*{sample}_*_L00*_R*_screen.txt"))

        for screen_file in screen_files:
            with screen_file.open("r") as f:
                df = pd.read_csv(f, sep="\t", header=1, skipfooter=2)

                # df.columns = df.columns.str.strip()

                for _, row in df.iterrows():
                    genome = row["Genome"]
                    one_hit = row["%One_hit_one_genome"]
                    multiple_hits = row["%Multiple_hits_one_genome"]

                    if genome not in sample_screens[sample]:
                        sample_screens[sample][genome] = {
                            "one_hit": [],
                            "multiple_hits": [],
                        }

                    sample_screens[sample][genome]["one_hit"].append(one_hit)
                    sample_screens[sample][genome]["multiple_hits"].append(
                        multiple_hits
                    )

    # Calculate mean percentages
    for sample, genomes in sample_screens.items():
        for genome, values in genomes.items():
            sample_screens[sample][genome]["one_hit"] = sum(values["one_hit"]) / len(
                values["one_hit"]
            )
            sample_screens[sample][genome]["multiple_hits"] = sum(
                values["multiple_hits"]
            ) / len(values["multiple_hits"])
            sample_screens[sample][genome]["total_hits"] = (
                sample_screens[sample][genome]["one_hit"]
                + sample_screens[sample][genome]["multiple_hits"]
            )

    return sample_screens


def save_to_csv(output_file, sample_totals, sample_screens, sample_flagstats):
    output_file = Path(output_file)

    fieldnames = (
        ["Sample"]
        + list(sample_totals[next(iter(sample_totals))].keys())
        + list(sample_flagstats[next(iter(sample_flagstats))].keys())
    )

    all_genomes = set()
    for sample_data in sample_screens.values():
        all_genomes.update(sample_data.keys())

    for genome in all_genomes:
        fieldnames.append(f"Screen {genome} One Hit Percent")
        fieldnames.append(f"Screen {genome} Multiple Hit Percent")
        fieldnames.append(f"Screen {genome} Total Hit Percent")

    with output_file.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for sample in sample_totals.keys():
            row = {"Sample": sample}
            row.update(sample_totals.get(sample, {}))
            row.update(sample_flagstats.get(sample, {}))

            for genome in all_genomes:
                row[f"Screen {genome} One Hit Percent"] = (
                    sample_screens.get(sample, {}).get(genome, {}).get("one_hit", "N/A")
                )
                row[f"Screen {genome} Multiple Hit Percent"] = (
                    sample_screens.get(sample, {})
                    .get(genome, {})
                    .get("multiple_hits", "N/A")
                )
                row[f"Screen {genome} Total Hit Percent"] = (
                    sample_screens.get(sample, {})
                    .get(genome, {})
                    .get("total_hits", "N/A")
                )

            writer.writerow(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract total reads from JSON files and screen data from TXT files, then save to CSV."
    )
    parser.add_argument(
        "--input_json_dir", help="Path to the directory containing JSON files."
    )
    parser.add_argument(
        "--input_screen_dir", help="Path to the directory containing screen text files."
    )
    parser.add_argument(
        "--input_flagstat_dir", help="Path to the directory containing flagstat files."
    )
    parser.add_argument("--output_csv", help="Path to the output CSV file.")

    args = parser.parse_args()

    sample_list = [
        # "23060404", "23091057", "23113261", "23113287", "24011925",
        # "ANT5530", "ANT5629", "ANT5639", "ANT5661", "ANT5670", "ANT5699", "ANT5797", "YOK3341"
        # "106264-001-112",
        # "106264-001-113",
        # "106264-001-114",
        # "106264-001-115",
        # "106264-001-116",
        # "106264-001-117",
        # "106264-001-118",
        # "106264-001-119",
        # "106264-001-120",
        # "106264-001-121",
        # "106264-001-122",
        # "106264-001-123",
        # "106264-001-124",
        # "106264-001-125",
        # "106264-001-126",
        # "106264-001-127",
        # "106264-001-128"
        "106264-002-077",
        "106264-002-078",
        "106264-002-079",
        "106264-002-080",
        "106264-002-081",
        "106264-002-082",
        "106264-002-083",
        "106264-002-084",
        "106264-002-085",
        "106264-002-086",
        "106264-002-087",
        "106264-002-088",
        "106264-002-089",
        "106264-002-090",
        "106264-002-091",
        "106264-002-092",
        "106264-002-093",
        "106264-002-094",
        "106264-002-095",
        "106264-002-096",
        "106264-002-097",
        "106264-002-098",
    ]

    results = extract_total_reads(args.input_json_dir, sample_list)
    screen_results = extract_screen_data(args.input_screen_dir, sample_list)
    flagstat_results = extract_flagstat_data(args.input_flagstat_dir, sample_list)
    save_to_csv(args.output_csv, results, screen_results, flagstat_results)
    print(f"CSV file '{args.output_csv}' has been created.")
