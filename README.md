# Antibody Repertoire Analysis Pipeline

This repository contains a standardized, Snakemake-based pipeline for analyzing antibody repertoire sequencing data. The pipeline processes raw FASTQ files, performs clonotype analysis, evaluates diversity metrics, and optionally compares immunized vs. unimmunized samples. It integrates MiXCR for clonotype assembly and custom Python scripts for downstream analysis.

## Features

- **MiXCR Processing**: Aligns and assembles raw FASTQ files to extract clonotypes for heavy (IGH) and kappa light (IGK) chains using MiXCR.
- **Diversity Analysis**: Computes clonality, Shannon index, Gini index, and CDR3 length distribution for each sample.
- **VDJ/VJ Combination Analysis**: Identifies top VDJ (heavy chain) and VJ (light chain) combinations and visualizes their usage.
- **Replicate and Group Comparison**: Analyzes replicate similarity and compares groups (e.g., immunized vs. unimmunized) using Jaccard index and statistical tests.
- **Snakemake Workflow**: Modular, reproducible pipeline with dynamic sample handling and configuration via `config.yaml`.
- **Visualization**: Generates plots for CDR3 length, gene usage, VDJ combinations, and clone expansion profiles.

## Prerequisites

- **Software**:
  - [Snakemake](https://snakemake.readthedocs.io/) (>=6.0.0)
  - [MiXCR](https://milaboratory.com/software/mixcr/) (ensure executable is in PATH)
  - Python (>=3.8)
- **Python Dependencies**:
  ```bash
  pip install pandas numpy scipy seaborn matplotlib scikit-bio matplotlib-venn
  ```
- **System Requirements**:
  - High-memory system (MiXCR requires ~500GB RAM for large datasets)
  - Linux/Unix environment (tested on Ubuntu)

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/biomaticss/Antibody-Discovery-Pipeline.git
   cd antibody-repertoire-analysis
   ```

2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
   Create a `requirements.txt` with:
   ```text
   pandas
   numpy
   scipy
   seaborn
   matplotlib
   scikit-bio
   matplotlib-venn
   ```

3. **Install MiXCR**:
   Follow the [MiXCR installation guide](https://milaboratory.com/software/mixcr/installation/) and ensure `mixcr` is in your PATH.

4. **Configure the Pipeline**:
   Edit `config.yaml` to specify input/output directories, sample names, and parameters (e.g., read count thresholds).

## Directory Structure

```
antibody-repertoire-analysis/
├── Snakefile                # Snakemake workflow
├── config.yaml              # Configuration file
├── MIXCR.sh                 # MiXCR processing script (optional; inline in Snakefile)
├── ana_modified.py          # Diversity and clonotype analysis
├── vdj_modified.py          # Immunized vs. unimmunized comparison
├── requirements.txt         # Python dependencies
└── README.md                # This file
```

## Usage

1. **Prepare Input Data**:
   - Place paired-end FASTQ files (`SAMPLEID_1.fq.gz`, `SAMPLEID_2.fq.gz`) in the input directory specified in `config.yaml`.
   - Update `config.yaml` with sample IDs and paths.

2. **Run the Pipeline**:
   - **Basic Analysis** (MiXCR + diversity):
     ```bash
     snakemake --cores 4 all
     ```
   - **With VDJ Comparison** (if immunized/unimmunized summaries exist):
     ```bash
     snakemake --cores 4 all_with_vdj
     ```

3. **Visualize Workflow**:
   ```bash
   snakemake --dag | dot -Tpng > dag.png
   ```

4. **Outputs**:
   - **MiXCR**: Clonotype TSV files (`SAMPLEID_VH_assemble_Clones_IGH.tsv`, `SAMPLEID_VK_assemble_Clones_IGK.tsv`) in `mixcr_output_dir`.
   - **Analysis**: CSV files, plots (CDR3 length, gene usage, VDJ combinations), and summary Excel (`top_combinations_summary.xlsx`) in `ana_output_dir`.
   - **VDJ Comparison**: CSV files (`top_vh_differences.csv`, `top_vk_differences.csv`) in `vdj_output_dir`.

5. **Clean Up**:
   ```bash
   snakemake clean
   ```

## Configuration (`config.yaml`)

Key fields:
- `input_dir`: Directory with raw FASTQ files.
- `mixcr_output_dir`: MiXCR output directory.
- `ana_output_dir`: Analysis output directory.
- `vdj_output_dir`: VDJ comparison output directory.
- `samples`: Dictionary of sample IDs and their MiXCR data directories.
- `unimmunized_summary` / `immunized_summary`: Excel files from prior runs for VDJ comparison.
- Parameters: `mixcr_memory`, `ana_readcount_threshold`, `ana_top_n`, `vdj_top_n`, etc.

Example:
```yaml
input_dir: "/path/to/rawdata"
mixcr_output_dir: "/path/to/MixCR"
ana_output_dir: "/path/to/Analysis"
vdj_output_dir: "/path/to/VDJ_Comparison"
samples:
  1175-20: "/path/to/MixCR"
  1175-25: "/path/to/MixCR"
mixcr_memory: "500g"
ana_readcount_threshold: 3
ana_top_n: 100
```

## Notes

- **MiXCR Customization**: If your data includes IGL (λ light chain), modify `Snakefile` to export IGL clonotypes using `mixcr exportClones --chains IGL`.
- **Scalability**: Use a cluster (e.g., SLURM) with `snakemake --cluster "sbatch"` for large datasets.
- **Debugging**: Run `snakemake --dry-run` to verify dependencies. Logs are saved in `{output_dir}/logs/`.

## Contributing

Contributions are welcome! Please:
1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/YourFeature`).
3. Commit changes (`git commit -m "Add YourFeature"`).
4. Push to the branch (`git push origin feature/YourFeature`).
5. Open a pull request.

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/biomaticss/Antibody-Discovery-Pipeline/blob/Anti-library/LICENSE) file for details.

## Contact

For issues or questions, open a GitHub issue or contact [hyut_2000@163.com](hyut_2000@163.com).
