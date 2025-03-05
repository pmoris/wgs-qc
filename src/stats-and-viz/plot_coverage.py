import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# import numpy as np
import argparse
import subprocess
from pathlib import Path
import warnings

matplotlib.use("Agg")  # Non-GUI backend (for saving images)

# TODO: depth < 1 needs to be set to 0 or 1 to avoid issues in log scale OR set max y-axis value in absolute values
# TODO: set limit to regions for species with too many regions, or concat all non-chromosome regions
# TODO: allow script to take in a samplesheet with filepaths as input (or with sample names to search for?)


def run_samtools_depth(bam_file, depth_file):
    """Runs `samtools depth` and saves output to a file."""
    print(f"Generating depth file for {bam_file} -> {depth_file}")

    result = subprocess.run(
        ["samtools", "depth", "-a", str(bam_file)], capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"Error running samtools depth: {result.stderr}")

    # with open(depth_file, "w") as f:
    with depth_file.open("w") as f:
        f.write(result.stdout)

    data = [line.split("\t") for line in result.stdout.strip().split("\n")]
    return pd.DataFrame(data, columns=["chrom", "pos", "depth"]).astype(
        {"pos": int, "depth": int}
    )


# def plot_depth(df, bin_size, facet, log_scale, output_file=None, width=15):
def plot_depth(df, output_file, width, **kwargs):
    """Plots the coverage depth."""

    # Compute chromosome offsets
    chrom_lengths = df.groupby("chrom")["pos"].max().to_dict()
    cumulative_pos, offset = {}, 0
    for chrom, length in chrom_lengths.items():
        cumulative_pos[chrom] = offset
        offset += length
    df["global_pos"] = df["pos"] + df["chrom"].map(cumulative_pos)

    bin_size = kwargs.get("bin_size")

    # Plot
    if kwargs.get("facet"):
        # Bin coverage
        df["bin"] = (df["pos"] // bin_size) * bin_size
        binned_df = df.groupby(["bin", "chrom"], as_index=False)["depth"].mean()

        unique_chroms = binned_df["chrom"].unique()

        if kwargs.get("log_scale"):
            binned_df["depth"] = binned_df["depth"] + 1

        fig, axes = plt.subplots(
            len(unique_chroms), 1, figsize=(width, 3 * len(unique_chroms)), sharex=False
        )

        if len(unique_chroms) == 1:
            axes = [axes]

        for ax, (chrom, chrom_df) in zip(axes, binned_df.groupby("chrom")):
            ax.plot(chrom_df["bin"], chrom_df["depth"], color="blue", linewidth=0.7)
            ax.set_xlim(
                [0, chrom_df["bin"].max()]
            )  # Set X-axis to range from 0 to max position for this chromosome
            ax.set_title(f"Coverage for {chrom}")
            ax.set_ylabel("Depth")
            ax.grid(True, linestyle="--", alpha=0.5)

            if kwargs.get("log_scale"):
                ax.set_yscale("log")
                ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

        plt.xlabel("Genomic Position")

    else:
        # Bin coverage
        df["bin"] = (df["global_pos"] // bin_size) * bin_size
        binned_df = df.groupby(["bin", "chrom"], as_index=False)["depth"].mean()

        if kwargs.get("log_scale"):
            binned_df["depth"] = binned_df["depth"] + 1

        plt.figure(figsize=(width, 5))
        plt.plot(binned_df["bin"], binned_df["depth"], color="blue", linewidth=0.7)

        # Add chromosome tick labels & alternating shading
        chrom_midpoints, prev_end = {}, 0
        for i, (chrom, start) in enumerate(cumulative_pos.items()):
            chrom_mid = start + (chrom_lengths[chrom] // 2)
            chrom_midpoints[chrom] = chrom_mid
            plt.axvspan(
                prev_end,
                start + chrom_lengths[chrom],
                color="gray",
                alpha=0.2 if i % 2 == 0 else 0,
            )
            prev_end = start + chrom_lengths[chrom]

        plt.xticks(
            list(chrom_midpoints.values()),
            labels=list(chrom_midpoints.keys()),
            rotation=45,
            ha="right",
        )
        plt.xlabel("Chromosome")
        plt.ylabel("Depth")
        plt.title("Coverage Plot Across Genome")
        plt.grid(True, linestyle="--", alpha=0.5)

        if kwargs.get("log_scale"):
            plt.yscale("log")
            plt.gca().yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            plt.ylabel("Depth (log10)")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def plot_cumulative_coverage(df, output_file=None, width=15, log_scale=False):
    """Plots the cumulative coverage distribution."""
    coverage_counts = df["depth"].value_counts().sort_index()
    cumulative_fraction = 1 - coverage_counts.sort_index().cumsum() / len(df)

    plt.figure(figsize=(width, 5))
    plt.plot(
        cumulative_fraction.index, cumulative_fraction, color="blue", linewidth=0.7
    )
    plt.xlabel("Coverage Depth")
    plt.ylabel("Fraction of Genome")
    plt.title("Cumulative Coverage Distribution")
    plt.grid(True, linestyle="--", alpha=0.5)

    if log_scale:
        plt.xscale("log")
        plt.gca().xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        plt.xlabel("Coverage Depth (log10)")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
    # TODO: write logic for when no output file is provided
    # else:
    #     plt.show()


def is_bam_file(input_file):
    "Check if input file is bam or binary"
    # if input_file.suffix == "bam" or input_file.suffix == "BAM":
    #     return True
    # else:
    #     try:
    #         with open(input_file, 'rb') as f:
    #             chunk = f.read(1024)
    #         return b'\x00' in chunk  # Presence of null bytes indicates binary
    #     except Exception as e:
    #         print(f"Error reading file: {e}")
    #         return False  # Assume non-binary in case of an error

    try:
        with open(input_file, "rb") as f:
            chunk = f.read(1024)
        chunk.decode("utf-8")  # Try decoding as UTF-8
        return False  # If decoding succeeds, it's a text file
    except UnicodeDecodeError:
        return True  # If decoding fails, it's a binary file


def plot_samtools_depth(
    output_file,
    bam_file=None,
    depth_file=None,
    bin_size=1000,
    width=15,
    facet=False,
    log_scale=False,
    chrom_filter=None,
    cumulative_coverage=False,
    overwrite=False,
    overwrite_depth=False,
):
    """
    Plots coverage depth from samtools depth output.

    Parameters:
        output_file (str): Output file locatoin to save the plot to.
        bam_file (str): Path to BAM file (optional, triggers samtools depth).
        depth_file (str): Path to the samtools depth output file or output location when generating from bam file (optional).
        bin_size (int): Size of bins for averaging coverage.
        width (int): Width of the plot.
        facet (bool): Whether to facet into subplots per chromosome.
        log_scale (bool): If True, uses log-scale for Y-axis.
        chrom_filter (str): If provided, only plot this chromosome.
        cumulative_coverage (bool): Whether to plot cumulative coverage distribution.
    """
    # Create or load depth data depending on supplied files
    if bam_file:
        # if not bam_file.exists():
        #     raise FileNotFoundError(f"BAM file not found: {bam_file}")
        if depth_file.exists() and not overwrite_depth:
            print(
                f"Skipping samtools depth since {depth_file} already exists. Use --overwrite-depth to force re-generation."
            )
            df = pd.read_csv(
                depth_file, sep="\t", header=None, names=["chrom", "pos", "depth"]
            )
        else:
            df = run_samtools_depth(bam_file, depth_file)
    elif depth_file:
        # if not depth_file.exists():
        #     raise FileNotFoundError(f"Depth file not found: {depth_file}")
        df = pd.read_csv(
            depth_file, sep="\t", header=None, names=["chrom", "pos", "depth"]
        )
    else:
        raise ValueError(
            "Either --input-depth or --input-bam must be provided (or --input-directory)."
        )

    if output_file.exists() and not overwrite:
        print(
            f"Skipping {output_file}, already exists. Use --overwrite to force re-generation."
        )
        return

    # Filter by chromosome if specified
    if chrom_filter:
        df = df[df["chrom"] == chrom_filter]

    if cumulative_coverage:
        plot_cumulative_coverage(df, output_file, width, log_scale)
    else:
        plot_depth(
            df=df,
            output_file=output_file,
            width=width,
            bin_size=bin_size,
            facet=facet,
            log_scale=log_scale,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Plot depth or cumulative coverage from BAM files or samtools depth output."
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file for an individual plot or path to directory when processing a directory of files.",
        required=True,
    )
    parser.add_argument(
        "--output-suffix",
        type=str,
        default=None,
        help="File suffix for output, defaults to depth or cumulative-coverage.",
    )
    parser.add_argument(
        "--input-dir", type=str, help="Input directory containing BAM or depth files."
    )
    parser.add_argument(
        "--input-pattern",
        type=str,
        # default="*.bam",
        help="Pattern to match files.",
    )
    parser.add_argument(
        "--input-bam", type=str, help="Path to a single input BAM file."
    )
    parser.add_argument(
        "--input-depth", type=str, help="Path to a single input depth file."
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing figures."
    )
    parser.add_argument(
        "--overwrite-depth", action="store_true", help="Overwrite existing depth files."
    )
    parser.add_argument(
        "--bin-size", type=int, default=1000, help="Bin size for averaging coverage."
    )
    parser.add_argument(
        "--width", type=int, default=15, help="Plot width (default: 15)."
    )
    parser.add_argument(
        "--facet", action="store_true", help="Create facetted subplots per chromosome."
    )
    parser.add_argument("--log", action="store_true", help="Use log-scale for depth.")
    parser.add_argument("--region", type=str, help="Select a specific region to plot.")
    parser.add_argument(
        "--cumulative-coverage",
        action="store_true",
        help="Plot cumulative coverage instead of per-region depth.",
    )
    args = parser.parse_args()

    if args.facet and args.cumulative_coverage:
        raise ValueError("Facet and cumulative coverage modes are not compatible.")

    suffix = args.output_suffix or (
        "cumulative_coverage" if args.cumulative_coverage else "depth"
    )

    if args.input_dir:
        input_dir = Path(args.input_dir)
        if not input_dir.exists():
            raise FileNotFoundError(f"Input directory not found: {input_dir}")
        if not args.input_pattern:
            print(
                "No input file pattern provided, defaulting to files matching '*.bam'."
            )

        if args.input_bam or args.input_depth:
            raise UserWarning(
                "Input bam or depth file should not be given when processing a directory of files."
            )

        pattern = args.input_pattern or ".bam"
        input_files = list(input_dir.glob("*" + pattern))

        if not input_files:
            raise FileNotFoundError(
                f"No input files found matching {pattern} in directory {input_dir}."
            )

        print(f"Processing all files in {input_dir} that match the pattern {pattern}.")

        output_dir = Path(args.output)
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)

        for file in input_files:
            output_file = (
                (output_dir / f"{file.stem}_{suffix}.png")
                if output_dir
                else file.with_suffix(f"_{suffix}.png")
            )

            if output_file.exists() and not args.overwrite:
                print(
                    f"Skipping creation of {output_file}, already exists. Use --overwrite to force re-generation."
                )
                continue

            if is_bam_file(file):
                # TODO: move testing logic to inside function?
                bam_file = Path(file)
                depth_file = Path(output_dir) / f"{file.stem}.depth"
                if not bam_file.exists():
                    raise FileNotFoundError(f"BAM file not found: {bam_file}")

            else:
                depth_file = Path(file)
                bam_file = None
                if not depth_file.exists():
                    raise FileNotFoundError(f"Depth file not found: {depth_file}")
                print(f"Using existing depth file: {depth_file}")

            print(f"Processing {file} → {output_file}")
            plot_samtools_depth(
                output_file=output_file,
                depth_file=depth_file,
                bam_file=bam_file,
                bin_size=args.bin_size,
                width=args.width,
                facet=args.facet,
                log_scale=args.log,
                chrom_filter=args.region,
                cumulative_coverage=args.cumulative_coverage,
                overwrite=args.overwrite,
                overwrite_depth=args.overwrite_depth,
            )
    else:
        if args.input_pattern:
            raise UserWarning(
                "Pattern should not be supplied when processing a single file."
            )

        if args.input_depth and args.input_bam:
            raise UserWarning(
                "Both an input BAM and depth file were provided, while only one should be used at a time."
            )

        output_file = Path(args.output)

        if output_file.exists() and not args.overwrite:
            print(
                f"Skipping creation of {output_file}, already exists. Use --overwrite to force re-generation."
            )
            return

        if not output_file.parent.exists():
            raise ValueError(
                f"Output file directory {output_file.parent} does not exist."
            )

        bam_file = None
        depth_file = None

        if args.input_bam:
            bam_file = Path(args.input_bam)
            if not bam_file.exists():
                raise FileNotFoundError(f"Input BAM file not found: {bam_file}")
            depth_file = output_file.parent / f"{bam_file.stem}.depth"

        if args.input_depth:
            depth_file = Path(args.input_depth)
            if not depth_file.exists():
                raise FileNotFoundError(f"Input depth file not found: {depth_file}")

        plot_samtools_depth(
            output_file=output_file,
            bam_file=bam_file,
            depth_file=depth_file,
            bin_size=args.bin_size,
            width=args.width,
            facet=args.facet,
            log_scale=args.log,
            chrom_filter=args.region,
            cumulative_coverage=args.cumulative_coverage,
            overwrite=args.overwrite,
            overwrite_depth=args.overwrite_depth,
        )


if __name__ == "__main__":
    main()


# TODO: switch back to function taking both bam and depth file (default None) and then choosing logic to follow based on which one is provided
# for batch mode, two different function calls will be necessary
