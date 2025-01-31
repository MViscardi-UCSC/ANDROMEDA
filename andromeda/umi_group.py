"""
umi_group.py
Marcus Viscardi,    January 31, 2025

This is really just a wrapper around the UMI_tools group command, but it should make it easier to run the command
and then plot the distribution of UMI group sizes.

# TODO: It would be nice to have some sort of inference regarding what edit distance threshold
        to use depending on the type and length of the UMI sequences being used.
        But that's a problem for another day. And if we have this script be easy enough to use
        I could write a wrapper to just try a bunch of different edit distances and plot the results
        to help the user decide what to use!!
"""
import argparse
import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm.auto import tqdm


def run_umi_tools_group(
        input_bam: Path,
        output_dir: Path,
        umi_tag: str = "uM",
        cell_tag: str = None,  # This can also be thought of a 2nd UMI! It would work well for RT tag and PCR tags!
        edit_dist: int = 1,
        per_contig: bool = True,  # This will merge all reads with matching contigs and UMIs into the same group!
        per_gene: bool = False,  # This will merge all reads with matching genes and UMIs into the same group!
        per_cell: bool = False  # This will require above selection and a matching cell barcode (or 2nd UMI) to group!!
) -> Path:
    """
    Runs `umi_tools group` on a BAM file with a single edit distance threshold.

    Args:
        input_bam (Path): Path to input BAM file (must have UMIs tagged).
        output_dir (Path): Directory to save grouped BAM and summary files.
        umi_tag (str): BAM tag where UMIs are stored (default: "u1").
        cell_tag (str): BAM tag for cell barcodes (if applicable). This can be used as a 2nd UMI! Good for RTvPCR!
        edit_dist (int): Edit distance threshold for UMI grouping.
        per_contig (bool): Merge all reads with matching contigs and UMIs into the same group
        per_gene (bool): Merge all reads with matching genes and UMIs into the same group
        per_cell (bool): Merge all reads with matching contigs/genes, cell barcodes and UMIs into the same group

    Returns:
        Path: Path to the grouped BAM file.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"📍 Running UMI grouping with edit distance {edit_dist}")

    grouped_bam = output_dir / input_bam.with_suffix(f".grouped_{edit_dist}dist.bam").name
    group_out_tsv = output_dir / input_bam.with_suffix(f".grouped_{edit_dist}dist.tsv").name
    group_log = output_dir / input_bam.with_suffix(f".grouped_{edit_dist}dist.log").name

    # Construct `umi_tools group` command
    group_call = [
        "umi_tools", "group",
        "-I", str(input_bam),
        "--group-out", str(group_out_tsv),
        "--method", "directional",  # This seemed to work well (and is default) but we could try others!
        "--edit-distance-threshold", str(edit_dist),
        "--log", str(group_log),
        "--extract-umi-method", "tag",
        "--umi-tag", umi_tag,
        "--output-bam", "--stdout", str(grouped_bam)
    ]

    # Optional parameters
    if per_contig:
        # group_call.append("--per-contig")  # Main thing we want to do for most of our custom/targeted sequencing!
        group_call.extend(["--per-gene", "--per-contig"])
    elif per_gene:
        group_call.append("--per-gene")
    if per_cell and cell_tag:
        group_call.extend(["--per-cell", "--cell-tag", cell_tag])
    elif per_cell or cell_tag:
        raise ValueError("per_cell requires a cell_tag to be provided! and vice versa!")

    # Run UMI grouping
    process = subprocess.Popen(group_call, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # Print stdout and stderr live
    for stdout_line in iter(process.stdout.readline, ""):
        print(stdout_line, end="")
    for stderr_line in iter(process.stderr.readline, ""):
        print(stderr_line, end="")
    
    process.stdout.close()
    process.stderr.close()
    process.wait()
    
    # subprocess.run(group_call check=True)  # TODO: Put this less crazy version back if it works better

    # Index BAM
    subprocess.run(["samtools", "index", str(grouped_bam)], check=True)

    print(f"✅ Grouped BAM saved: {grouped_bam}")
    return grouped_bam


def plot_umi_distribution(grouping_tsv: Path, output_dir: Path):
    """
    Plots the distribution of UMI group sizes from `umi_tools group` output.

    Args:
        grouping_tsv (Path): TSV file containing grouped UMIs.
        output_dir (Path): Directory to save plots.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    edit_dist = grouping_tsv.stem.split("_")[-1].replace("dist", "")

    df = pd.read_csv(grouping_tsv, sep="\t")
    df = df[["final_umi", "final_umi_count", "unique_id"]].drop_duplicates()

    plt.figure(figsize=(6, 4))
    sns.histplot(df["final_umi_count"], bins=30, log_scale=True, kde=True, color="blue")
    plt.title(f"UMI Group Sizes (Edit Distance: {edit_dist})")
    plt.xlabel("UMI Count")
    plt.ylabel("Frequency")
    plt.grid(True)

    output_plot = output_dir / f"umi_distribution_{edit_dist}.png"
    plt.savefig(output_plot, dpi=300)
    print(f"📊 Saved plot: {output_plot}")
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Group UMIs using `umi_tools group` with a single edit distance.")

    parser.add_argument("input_bam", type=Path,
                        help="Path to input BAM file (must have UMIs tagged).")
    parser.add_argument("output_dir", type=Path,
                        help="Output directory to save grouped BAM and plots.")
    parser.add_argument("--umi-tag", type=str, default="uM",
                        help="BAM tag for UMIs (default: 'uM').")
    parser.add_argument("--cell-tag", type=str, default=None,
                        help="BAM tag for 'cell' barcodes (these can also be considered a second UMI ).")
    parser.add_argument("--edit-dist", type=int, required=True,
                        help="Edit distance threshold for UMI grouping.")
    parser.add_argument("--per-contig", action="store_true",
                        help="Group per contig. (highly recommended if your data is "
                             "targeted sequencing or exogenous RNA)")
    parser.add_argument("--per-gene", action="store_true",
                        help="Group per gene.")
    parser.add_argument("--per-cell", action="store_true",
                        help="Group per cell barcode (requires --cell-tag).")
    parser.add_argument("--plot", action="store_true",
                        help="Plot UMI group size distribution.")

    return parser


def group_umis(args):
    # Run UMI grouping
    grouped_bam = run_umi_tools_group(
        input_bam=args.input_bam,
        output_dir=args.output_dir,
        umi_tag=args.umi_tag,
        cell_tag=args.cell_tag,
        edit_dist=args.edit_dist,
        per_contig=args.per_contig,
        per_gene=args.per_gene,
        per_cell=args.per_cell,
    )

    # Generate UMI grouping statistics
    tsv_file = grouped_bam.with_suffix(".tsv")
    if args.plot:  # TODO: Add a "just-plot" option to the CLI, so we can run this script on existing data!
        plot_umi_distribution(tsv_file, args.output_dir / "plots")

def main():
    parser = parse_args()
    args = parser.parse_args()

    # Run UMI grouping
    grouped_bam = run_umi_tools_group(
        input_bam=args.input_bam,
        output_dir=args.output_dir,
        umi_tag=args.umi_tag,
        cell_tag=args.cell_tag,
        edit_dist=args.edit_dist,
        per_contig=args.per_contig,
        per_gene=args.per_gene,
        per_cell=args.per_cell,
    )

    # Generate UMI grouping statistics
    tsv_file = grouped_bam.with_suffix(".tsv")
    plot_umi_distribution(tsv_file, args.output_dir / "plots")


if __name__ == "__main__":
    main()
