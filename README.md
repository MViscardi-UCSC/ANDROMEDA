```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃            _   _ _____  _____   ____  __  __ ______ _____            ┃
┃      /\   | \ | |  __ \|  __ \ / __ \|  \/  |  ____|  __ \   /\      ┃
┃     /  \  |  \| | |  | | |__) | |  | | \  / | |__  | |  | | /  \     ┃
┃    / /\ \ | . ` | |  | |  _  /| |  | | |\/| |  __| | |  | |/ /\ \    ┃
┃   / ____ \| |\  | |__| | | \ \| |__| | |  | | |____| |__| / ____ \   ┃
┃  /_/    \_\_| \_|_____/|_|  \_\\____/|_|  |_|______|_____/_/    \_\  ┃
┃                                                                      ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](https://www.python.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![GitHub License](https://img.shields.io/github/license/MViscardi-UCSC/ANDROMEDA)](https://github.com/MViscardi-UCSC/ANDROMEDA/blob/master/LICENSE)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/t/MViscardi-UCSC/ANDROMEDA)](https://github.com/MViscardi-UCSC/ANDROMEDA/commits/master/)
[![GitHub contributors](https://img.shields.io/github/contributors/MViscardi-UCSC/ANDROMEDA)](https://github.com/MViscardi-UCSC/ANDROMEDA/commits/master/)
[![GitHub Tag](https://img.shields.io/github/v/tag/MViscardi-UCSC/ANDROMEDA?logo=github)](https://github.com/MViscardi-UCSC/ANDROMEDA/tags)
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/MViscardi-UCSC/ANDROMEDA/python-package.yml?logo=github&label=CI&link=https%3A%2F%2Fgithub.com%2FMViscardi-UCSC%2FANDROMEDA%2Factions%2Fworkflows%2Fpython-package.yml)](https://github.com/MViscardi-UCSC/ANDROMEDA/actions/workflows/python-package.yml)


**Alignment-based Nucleotide Detection and Read Optimization for Mapping Errors, Deaminations, and Alterations**

***

## Overview
ANDROMEDA is a bioinformatics tool designed to process mapped sequencing reads and extract information from ambiguous nucleotide positions. The tool identifies barcode and UMI sequences from predefined mapping positions, assigns tags to BAM files, and detects mismatched bases resulting from RNA modifications (e.g., deamination by TadA or other base editors). It also enables UMI grouping, consensus calling, and error annotation to facilitate high-confidence sequence reconstruction.


## Key Features
- **Map-based Barcode & UMI Extraction** – Extracts nucleotide sequences from user-defined ambiguous positions in mapped reads.
- **BAM Tagging** – Annotates BAM files with extracted barcode/UMI sequences for downstream analysis.
- **Modification Detection** – Identifies mismatched bases caused by base modifications (e.g., deamination).
- **UMI Grouping & Consensus Calling** – Leverages [UMI-tools](https://github.com/CGATOxford/UMI-tools) to cluster reads and generate consensus sequences.
- **Error & Confidence Tagging** – Adds metadata on sequence confidence, ambiguous bases, and read support.
- **Modular Design** – Individual processing steps can be run independently, allowing seamless integration into other workflows.

## Running ANDROMEDA with uv
To run ANDROMEDA with uv, you will need to have the following installed:
- **uv** ([GitHub link](https://github.com/astral-sh/uv))
- **samtools** ([GitHub link](https://github.com/samtools/samtools))

The `uvx` command will do the rest (gathering the correct python version and all other dependencies)! ***WOW!***

To run the full pipeline, you can use the following command:
```bash
uvx --from git+https://github.com/MViscardi-UCSC/ANDROMEDA andromeda run-all --help
```

## Full Installation
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
3. consensus.py will take in the output from umi_group.py and generate consensus sequences for each UMI group. It will output a new BAM file of the consensus sequences.

More specific instructions for each module can be found below, with extended details available with the `--help` flag.

### 0. Pick Reference Positions for UMIs
```bash
python andromeda ref_pos_picker <ref.fasta> <output_parent_directory>
```
### 1. Extract Barcodes and UMIs from Mapped Reads
```bash
python andromeda extract <ref.fasta> <mapped.bam> <output_parent_directory>
```

### 2. Group UMIs Using UMI-Tools
```bash
python andromeda umi_group <tagged.bam (from extract step)> <output_parent_directory>
```

### 3. Create Consensus Sequences
```bash
python andromeda consensus <ref.fasta> <grouped.bam (from group step)> <output_parent_directory>
```

### Or Run all the steps together!:
```bash
python andromeda run-all <ref.fasta> <mapped.bam> <output_parent_directory>
```



## Dependencies
- `Python 3.11+`
- `samtools`
- `UMI-tools`
- `pandas`
- `biopython`
- `seaborn` and `matplotlib` for plotting *(~~optional~~)* ***TODO: Make this optional***


## Future Features
- Support for base modifications, not just mismatches due to deamination
- Integration with alternative consensus-calling algorithms
- Enhanced confidence scoring for ambiguous base calls

## License
MIT License. See `LICENSE` for details.

## Contributors
Developed by Marcus Viscardi and Liam Tran in the Arribere Lab at UCSC. Contributions welcome!

## Contact
For questions or feedback, please contact marcus.viscardi@gmail.com or open an issue on GitHub.
