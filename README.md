# ANDROMEDA
**Alignment-based Nucleotide Detection and Read Optimization for Mapping Errors, Deaminations, and Alterations**

<br>
<br>
<br>

Within the andromeda folder I have the first couple modules planned out:
1. extract.py will take in a mapped BAM file, a reference genome, and a list of positions to extract from the reads. It will output a new BAM file with the extracted sequences in the read tags.
2. umi_group.py will take in the output from extract.py and group reads by UMI. It will output a new BAM file with the UMI sequences in the read tags. (This will leverage UMI-Tools)
3. consensus.py will take in the output from umi_group.py and generate consensus sequences for each UMI group. It will output a new BAM file with the consensus sequences in the read tags.

*ChatGPT ate with this one... Nothing below this line is implimented yet, but it's a nice goal to aspire to!!*
***


## Overview
ANDROMEDA is a modular bioinformatics tool designed to process mapped sequencing reads and extract information from ambiguous nucleotide positions. The tool identifies barcode and UMI sequences from predefined mapping positions, assigns tags to BAM files, and detects mismatched bases resulting from RNA modifications (e.g., deamination by ADAR or other base editors). It also enables UMI grouping, consensus calling, and error annotation to facilitate high-confidence sequence reconstruction.


## Key Features
- **Map-based Barcode & UMI Extraction** – Extracts nucleotide sequences from user-defined ambiguous positions in mapped reads.
- **BAM Tagging** – Annotates BAM files with extracted barcode/UMI sequences for downstream analysis.
- **Modification Detection** – Identifies mismatched bases caused by base modifications (e.g., deamination, chemical treatments).
- **UMI Grouping & Consensus Calling** – Leverages UMI-tools to cluster reads and generate consensus sequences.
- **Error & Confidence Tagging** – Adds metadata on sequence confidence, ambiguous bases, and read support.
- **Modular Design** – Individual processing steps can be run independently, allowing seamless integration into other workflows.

## Installation
ANDROMEDA requires Python 3.8+ and dependencies such as `pysam`, `UMI-tools`, and `samtools`. Install using:

```bash
# Clone the repository
git clone https://github.com/your-repo/ANDROMEDA.git
cd ANDROMEDA

# Create a virtual environment (optional but recommended)
python -m venv andromeda_env
source andromeda_env/bin/activate  # On Windows use: andromeda_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage
ANDROMEDA consists of multiple modules that can be run independently or as part of a pipeline.

### 1. Extract Barcodes and UMIs from Mapped Reads
```bash
python andromeda_extract.py --input mapped_reads.bam --ref reference.fasta --output tagged_reads.bam
```

### 2. Group UMIs and Generate Consensus Sequences
```bash
python andromeda_consensus.py --input tagged_reads.bam --output consensus_reads.bam
```

### 3. Identify Mismatched Bases (RNA Modifications)
```bash
python andromeda_mismatch.py --input consensus_reads.bam --output modified_reads.bam
```

## Example Workflow
To run the full pipeline from mapped reads to consensus sequences:
```bash
python andromeda_extract.py --input mapped_reads.bam --ref reference.fasta --output tagged_reads.bam
python andromeda_consensus.py --input tagged_reads.bam --output consensus_reads.bam
python andromeda_mismatch.py --input consensus_reads.bam --output modified_reads.bam
```

## Output
- **Tagged BAM File** – Includes extracted barcodes/UMIs in specified BAM tags.
- **Consensus BAM File** – Contains grouped and collapsed UMI reads with high-confidence sequences.
- **Modification-Tagged BAM File** – Annotates potential RNA modification sites and ambiguous base calls.

## Dependencies
- `Python 3.8+`
- `pysam`
- `UMI-tools`
- `samtools`
- `numpy`
- `pandas`

## Future Features
- Support for additional base modifications
- Integration with alternative consensus-calling algorithms
- Enhanced confidence scoring for ambiguous base calls

## License
MIT License. See `LICENSE` for details.

## Contributors
Developed by **[Your Name / Lab]**. Contributions welcome!

## Contact
For questions or feedback, please contact **[your email]** or open an issue on GitHub.
