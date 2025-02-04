# ANDROMEDA
**Alignment-based Nucleotide Detection and Read Optimization for Mapping Errors, Deaminations, and Alterations**

***

## Overview
ANDROMEDA is a modular bioinformatics tool designed to process mapped sequencing reads and extract information from ambiguous nucleotide positions. The tool identifies barcode and UMI sequences from predefined mapping positions, assigns tags to BAM files, and detects mismatched bases resulting from RNA modifications (e.g., deamination by ADAR or other base editors). It also enables UMI grouping, consensus calling, and error annotation to facilitate high-confidence sequence reconstruction.


## Key Features
- **Map-based Barcode & UMI Extraction** – Extracts nucleotide sequences from user-defined ambiguous positions in mapped reads.
- **BAM Tagging** – Annotates BAM files with extracted barcode/UMI sequences for downstream analysis.
- **Modification Detection** – Identifies mismatched bases caused by base modifications (e.g., deamination).
- **UMI Grouping & Consensus Calling** – Leverages UMI-tools to cluster reads and generate consensus sequences.
- **Error & Confidence Tagging** – Adds metadata on sequence confidence, ambiguous bases, and read support.
- **Modular Design** – Individual processing steps can be run independently, allowing seamless integration into other workflows.

## Installation
ANDROMEDA requires Python 3.8+ and dependencies such as `pysam`, `UMI-tools`, and `samtools`. Install using:

```bash
# Clone the repository
git clone https://github.com/MViscardi-UCSC/ANDROMEDA.git
cd ANDROMEDA

# Create a virtual environment (optional but recommended)
python -m venv andromeda_env
source andromeda_env/bin/activate  # On Windows use: andromeda_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage
ANDROMEDA consists of multiple modules that can be run independently ~~or as part of a pipeline~~ (coming soon).

Within the andromeda folder we currently have these modules working:

0. ref_pos_picker.py will help users pick correct coordinates for their UMIs in their reference
1. extract.py will take in a mapped BAM file, a reference genome, and a list of positions to extract from the reads. It will output a new BAM file with the extracted sequences in the read tags.
2. umi_group.py will take in the output from extract.py and group reads by UMI. It will output a new BAM file with the UMI sequences in the read tags. (This will leverage UMI-Tools)
3. consensus.py will take in the output from umi_group.py and generate consensus sequences for each UMI group. It will output a new ~~BAM file with the consensus sequences in the read tags~~ TSV file.


### 0. Pick Reference Positions for UMIs
```bash
python -m andromeda ref-pos-picker <ref.fasta>
```

### 1. Extract Barcodes and UMIs from Mapped Reads
```bash
python -m andromeda extract <reads.bam> <ref.fasta> <umi_pos_from_ref-pos-picker.tsv>
```

### 2. Group UMIs Using UMI-Tools
```bash
python -m andromeda group-umis <umi_tagged_reads.bam> <output_dir>
```

### 3. Create Consensus Sequences
```bash
python -m andromeda consensus <umi_grouped_reads.bam> <ref.fasta> <output_dir>
```

## Dependencies
- `Python 3.8+`
- `pysam`
- `UMI-tools`
- `samtools`
- `numpy`
- `pandas`

## Future Features
- 
- Support for additional base modifications
- Integration with alternative consensus-calling algorithms
- Enhanced confidence scoring for ambiguous base calls

## License
MIT License. See `LICENSE` for details.

## Contributors
Developed by Marcus Viscardi and Liam Tran in the Arribere Lab at UCSC. Contributions welcome!

## Contact
For questions or feedback, please contact marcus.viscardi@gmail.com or open an issue on GitHub.
