# ANDROMEDA Development Plan  

## 📌 Overview  
ANDROMEDA (**Alignment-based Nucleotide Detection and Read Optimization for Mapping Errors, Deaminations, and Alterations**) is a bioinformatics tool for:  
- Extracting barcode & UMI sequences from mapped BAM files.  
- Grouping reads by UMI (using UMI-Tools).  
- Generating consensus sequences for each UMI group.  
- Identifying base modifications caused by deaminases or other treatments.  

This document outlines the development plan for implementing ANDROMEDA as a **modular, maintainable, and scalable** Python package.

---

## 🏗️ Phase 1: Project Setup & Structure  
- [X] **1. Initialize the project repository**  
   - Create the `ANDROMEDA/` folder and subdirectories:  
     ```
     ANDROMEDA/
     ├── andromeda/          # Main package
     │   ├── __init__.py
     │   ├── extract.py       # Barcode & UMI extraction
     │   ├── umi_group.py     # UMI clustering (UMI-Tools)
     │   ├── consensus.py     # Consensus sequence generation
     │   ├── io_utils.py      # File handling utilities
     │   ├── bam_utils.py     # BAM file processing
     │   ├── config.py        # Configuration management
     │   ├── logging_utils.py # Logging support
     │   ├── cli.py           # Command-line interface
     ├── scripts/             # Helper scripts
     ├── tests/               # Unit & integration tests
     ├── docs/                # Documentation
     ├── examples/            # Example input/output
     ├── requirements.txt     # Dependencies
     ├── setup.py             # Installation script
     ├── README.md            # Project overview
     ├── LICENSE              # License
     ├── .gitignore           # Ignore unnecessary files
     ```
   - Add `README.md` and `.gitignore`.  

- [X] **2. Set up a virtual environment (if not already done)**  
   - Activate existing virtual environment.  
   - Install required dependencies:  
     ```bash
     pip install pysam umi-tools
     ```

- [X] **3. Convert Jupyter Notebooks into Modular Python Scripts**  
   - Identify reusable functions from notebooks.  
   - Break them into **extract.py, umi_group.py, and consensus.py**.  
   - Replace hardcoded file paths with function arguments.  

---

## 🚀 Phase 2: Core Module Development  

- [X] **4. Implement `extract.py`** (Barcode & UMI Extraction)  
   - **Inputs**:  
     - Mapped **BAM file**  
     - Reference **genome FASTA**  
     - **Positions** list for barcode/UMI extraction  
   - **Process**:  
     - Extract sequences from specified positions.  
     - Tag extracted sequences into BAM read tags.  
     - Use **pysam** for BAM file handling.  
   - **Outputs**: New BAM file with extracted sequences in read tags.  

- [X] **5. Implement `umi_group.py`** (UMI Grouping)  
   - **Inputs**: BAM file with UMI tags (from `extract.py`).  
   - **Process**:  
     - Use **UMI-Tools** to cluster reads by UMI.  
     - Assign UMI groups and tag BAM reads.  
   - **Outputs**: BAM file with UMI group tags (`UG`).  

- [X] **6. Implement `consensus.py`** (Consensus Sequence Generation)  
   - **Inputs**: BAM file with UMI groups (from `umi_group.py`).  
   - **Process**:  
     - Collapse reads into consensus sequences per UMI group.  
     - Add **metadata tags** (member count, confidence scores).  
   - **Outputs**: BAM file with consensus sequences in read tags.  
   - **TODO**: Get the final output into BAM file format, currently we output consensus sequences in TSV!!

---

## ⚙️ Phase 3: CLI & Usability  

- [x] **7. Implement Command-Line Interface (`cli.py`)**  
   - Users should be able to run individual steps:  
     ```bash
     python -m andromeda.cli extract --input mapped.bam --ref genome.fa --positions pos.txt --output extracted.bam
     ```
   - Use `argparse` to manage arguments.  

- [ ] **8. Implement `config.py` for Config Management**  
   - Store default parameters in **JSON/YAML**.  
   - Allow user overrides via CLI options.  

- [ ] **9. Add Logging (`logging_utils.py`)**  
   - Implement logging to track progress/errors.  
   - Replace `print()` statements with structured logs.  

---

## 🧪 Phase 4: Testing & Documentation  

- [ ] **10. Implement Unit Tests (`tests/`)**  
   - Use **pytest** for function-level tests.  
   - Small test BAM/reference datasets in `data/`.  

- [ ] **11. Expand Documentation (`docs/`)**  
   - `usage.md`: How to use each module.  
   - `architecture.md`: Software design decisions.  

- [ ] **12. Package the Tool (`setup.py`)**  
   - Make it **pip-installable**:  
     ```bash
     pip install .
     ```

---

## 🔮 Future Enhancements (Post-Initial Release)  
- **Parallel Processing:** Speed up BAM handling.  
- **Advanced Filtering:** Quality-based UMI selection.  
- **Snakemake/Nextflow Workflow:** For large-scale datasets.  

---
