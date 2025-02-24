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
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm.auto import tqdm


def run_umi_tools_group(
        tagged_bam: Path,
        output_parent_dir: Path,
        umi_tag: str = "uM",
        cell_tag: str = None,  # This can also be thought of a 2nd UMI! It would work well for RT tag and PCR tags!
        edit_dist: int = -1,
        edit_frac: float = 0.0,
        per_contig: bool = True,  # This will merge all reads with matching contigs and UMIs into the same group!
        per_gene: bool = False,  # This will merge all reads with matching genes and UMIs into the same group!
        per_cell: bool = False  # This will require above selection and a matching cell barcode (or 2nd UMI) to group!!
) -> Path:
    """
    Runs `umi_tools group` on a BAM file with a single edit distance threshold.

    Args:
        tagged_bam (Path): Path to input BAM file (must have UMIs tagged).
        output_parent_dir (Path): Directory to save grouped BAM and summary files.
        umi_tag (str): BAM tag where UMIs are stored (default: "u1").
        cell_tag (str): BAM tag for cell barcodes (if applicable). This can be used as a 2nd UMI! Good for RTvPCR!
        edit_dist (int): Edit distance threshold for UMI grouping.
        edit_frac (float): Fractional edit distance threshold for UMI grouping.
        per_contig (bool): Merge all reads with matching contigs and UMIs into the same group
        per_gene (bool): Merge all reads with matching genes and UMIs into the same group
        per_cell (bool): Merge all reads with matching contigs/genes, cell barcodes and UMIs into the same group

    Returns:
        Path: Path to the grouped BAM file.
    """
    
    assert output_parent_dir.exists(), f"Output directory does not exist: {output_parent_dir}"
    output_dir = output_parent_dir / "grouped"
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    grouped_bam = output_dir / tagged_bam.with_suffix(f".grouped_{edit_dist}dist.bam").name
    group_out_tsv = output_dir / tagged_bam.with_suffix(f".grouped_{edit_dist}dist.tsv").name
    group_log = output_dir / tagged_bam.with_suffix(f".grouped_{edit_dist}dist.log").name
    
    # Let's calculate the edit distance fraction if it wasn't provided
    if edit_frac <= 0 and edit_dist <= 0:
        # if this is the case we'll just default to 0
        print("No edit distance or fraction provided, defaulting to 0 (perfect match only).")
        edit_dist = 0
    elif edit_frac > 0 > edit_dist:
        # In this case we need to grab an example UMI and calculate the edit distance
        with pysam.AlignmentFile(tagged_bam, "rb") as bam:
            for read in bam:
                umi = read.get_tag(umi_tag)
                if umi:
                    break
        edit_dist = int(len(umi) * edit_frac)
        print(f"Calculated edit distance: {edit_dist} from fraction: {edit_frac}"
              f" and UMI length: {len(umi)}")
    elif edit_frac > 0 and edit_dist > 0:
        print("Both edit distance and fraction provided, using edit distance.")

    print(f"📍 Running UMI grouping with edit distance {edit_dist} "
          f"(this will not produce any outputs until complete, please be patient!)")
    # Construct `umi_tools group` command
    group_call = [
        "umi_tools", "group",
        "-I", str(tagged_bam),
        "--group-out", str(group_out_tsv),
        "--method", "directional",  # This seemed to work well (and is default) but we could try others!
        "--edit-distance-threshold", str(edit_dist),
        "--log", str(group_log),
        "--extract-umi-method", "tag",
        "--umi-tag", umi_tag,
        "--output-bam", "--stdout", str(grouped_bam)
    ]

    # Optional parameters
    if per_gene and not per_contig:
        group_call.append("--per-gene")
    elif per_gene and per_contig:
        group_call.extend(["--per-gene", "--per-contig"])
    elif per_contig:
        print("🔍 --per-contig grouping option selected, adding --per-gene (umi-tools requires this).")
        group_call.extend(["--per-gene", "--per-contig"])
    else:
        print("🔍 No grouping option selected, defaulting to per-contig.")
        group_call.extend(["--per-gene", "--per-contig"])
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

    parser.add_argument("tagged_bam", type=Path,
                        help="Path to input BAM file (must have UMIs tagged).")
    parser.add_argument("output_parent_dir", type=Path,
                        help="Parent directory to make a new directory inside to save outputs.")
    parser.add_argument("--load-umi-tag", type=str, default="uM",
                        help="BAM tag to load UMIs from (default: 'uM').")
    parser.add_argument("--cell-tag", type=str, default=None,
                        help="BAM tag for 'cell' barcodes (these can also be considered a second UMI ).")
    parser.add_argument("--grouping-edit-dist", type=int, default=-1,
                        help="Edit distance threshold for UMI grouping "
                             "(defaults to 0 if this or frac not changed).")
    parser.add_argument("--grouping-edit-frac", type=float, default=-1,
                        help="Fractional edit distance threshold for UMI grouping "
                             "(defaults to 0 if this or dist not changed).")
    parser.add_argument("--per-contig", action="store_true",
                        help="Group per contig. (highly recommended if your data is "
                             "targeted sequencing or exogenous RNA, will be defaulted "
                             "to if no other option is selected)")
    parser.add_argument("--per-gene", action="store_true",
                        help="Group per gene.")
    parser.add_argument("--per-cell", action="store_true",
                        help="Group per cell barcode (requires --cell-tag).")
    parser.add_argument("--grouping-plot", action="store_true",
                        help="Plot UMI group size distribution.")

    return parser


def dependencies():
    return {
        "tagged_bam": "extract.tagged_bam",
        "output_parent_dir": "ref_pos_picker.output_parent_dir",
    }


def group_umis(args):
    # Run UMI grouping
    grouped_bam = run_umi_tools_group(
        tagged_bam=args.tagged_bam,
        output_parent_dir=args.output_parent_dir,
        umi_tag=args.load_umi_tag,
        cell_tag=args.cell_tag,
        edit_dist=args.grouping_edit_dist,
        edit_frac=args.grouping_edit_frac,
        per_contig=args.per_contig,
        per_gene=args.per_gene,
        per_cell=args.per_cell,
    )
    
    # Generate UMI grouping statistics
    tsv_file = grouped_bam.with_suffix(".tsv")
    if args.grouping_plot:  # TODO: Add a "just-plot" option to the CLI, so we can run this script on existing data!
        plot_umi_distribution(tsv_file, args.output_parent_dir / "grouping" / "plots")


def main():
    parser = parse_args()
    args = parser.parse_args()
    pipeline_main(args)


def pipeline_main(args):
    # Run UMI grouping
    grouped_bam = run_umi_tools_group(
        tagged_bam=args.tagged_bam,
        output_parent_dir=args.output_parent_dir,
        umi_tag=args.load_umi_tag,
        cell_tag=args.cell_tag,
        edit_dist=args.grouping_edit_dist,
        edit_frac=args.grouping_edit_frac,
        per_contig=args.per_contig,
        per_gene=args.per_gene,
        per_cell=args.per_cell,
    )

    # Generate UMI grouping statistics
    tsv_file = grouped_bam.with_suffix(".tsv")
    if args.grouping_plot:  # TODO: Add a "just-plot" option to the CLI, so we can run this script on existing data!
        plot_umi_distribution(tsv_file, args.output_parent_dir / "grouping" / "plots")
        # "grouped-bam": "umi_group.grouped_bam",
        # "ref-fasta": "Reference FASTA file",
        # "output-parent-dir": "Output directory for consensus sequences",
        # "--min-group-size": "Minimum reads per UMI group",
        # "--consensus-plot": "Plot consensus quality"
    pass_fwd_dict = {
        "grouped_bam": grouped_bam,
        "grouping_tsv": tsv_file,
        "output_parent_dir": args.output_parent_dir,
    }
    return pass_fwd_dict


if __name__ == "__main__":
    main()
