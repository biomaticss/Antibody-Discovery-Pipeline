# Snakefile for Antibody Repertoire Analysis Pipeline
import os
import glob
import yaml

configfile: "config.yaml"

# Load config
INPUT_DIR = config["input_dir"]
MIXCR_OUTPUT_DIR = config["mixcr_output_dir"]
ANA_OUTPUT_DIR = config["ana_output_dir"]
VDJ_OUTPUT_DIR = config["vdj_output_dir"]
SAMPLES = config["samples"]
UNIMMUNIZED_SUMMARY = config["unimmunized_summary"]
IMMUNIZED_SUMMARY = config["immunized_summary"]
UNIMM_VH_SHEETS = config["unimm_vh_sheets"]
UNIMM_VK_SHEETS = config["unimm_vk_sheets"]
IMM_VH_SHEETS = config["imm_vh_sheets"]
IMM_VK_SHEETS = config["imm_vk_sheets"]

# Create directories
os.makedirs(MIXCR_OUTPUT_DIR + "/logs", exist_ok=True)
os.makedirs(ANA_OUTPUT_DIR + "/logs", exist_ok=True)
os.makedirs(VDJ_OUTPUT_DIR + "/logs", exist_ok=True)

# Wildcard constraints
wildcard_constraints:
    sample = "|".join(SAMPLES.keys())

# Target rules
rule all:
    input:
        expand(ANA_OUTPUT_DIR + "/top_combinations_summary.xlsx"),
        expand(ANA_OUTPUT_DIR + "/{sample}_analysis_complete.txt", sample=SAMPLES.keys()),

rule all_with_vdj:
    input:
        rules.all.input,
        VDJ_OUTPUT_DIR + "/top_vh_differences.csv",
        VDJ_OUTPUT_DIR + "/top_vk_differences.csv",

# Rule 1: MixCR processing (parallel per sample)
rule mixcr_process:
    input:
        r1 = ancient(INPUT_DIR + "/{sample}_1.fq.gz"),
        r2 = ancient(INPUT_DIR + "/{sample}_2.fq.gz")
    output:
        tsv_igh = MIXCR_OUTPUT_DIR + "/{sample}_VH_assemble_Clones_IGH.tsv",
        tsv_igk = MIXCR_OUTPUT_DIR + "/{sample}_VK_assemble_Clones_IGK.tsv"
    log:
        MIXCR_OUTPUT_DIR + "/logs/{sample}_mixcr.log"
    resources:
        mem = config["mixcr_memory"]
    shell:
        """
        # Inline MixCR commands (adapted from MIXCR.sh for single sample)
        OUTPUT_VDJCA="{output.tsv_igh}.vdjca"
        OUTPUT_CLNS="{output.tsv_igh}.clns"
        
        mixcr -Xmx{resources.mem} align \\
            --preset generic-amplicon \\
            --species hsa \\
            --rna \\
            --assemble-clonotypes-by CDR3 \\
            --rigid-left-alignment-boundary \\
            --rigid-right-alignment-boundary C \\
            -OmergerParameters.minimalOverlap=10 \\
            -OmergerParameters.minimalIdentity=0.9 \\
            {input.r1} {input.r2} $OUTPUT_VDJCA 2>> {log}
        
        mixcr -Xmx{resources.mem} assemble \\
            --assemble-clonotypes-by "{{FR1Begin:FR4End}}" \\
            -OseparateByC=true \\
            -ObadQualityThreshold=0 \\
            $OUTPUT_VDJCA $OUTPUT_CLNS 2>> {log}
        
        mixcr -Xmx{resources.mem} exportClones \\
            $OUTPUT_CLNS \\
            {output.tsv_igh} 2>> {log}
        
        # For IGK (VK): Assume same process, but adjust if needed for light chain
        mixcr -Xmx{resources.mem} exportClonesLightChains \\
            $OUTPUT_CLNS \\
            {output.tsv_igk} 2>> {log}
        
        # Cleanup
        rm -f $OUTPUT_VDJCA $OUTPUT_CLNS
        
        # Check output
        if [[ ! -s {output.tsv_igh} || ! -s {output.tsv_igk} ]]; then
            echo "Error: Empty output for {wildcards.sample}" >> {log}
            exit 1
        fi
        echo "MixCR completed for {wildcards.sample}" >> {log}
        """

# Rule 2: Ana analysis per sample
rule ana_analyze_sample:
    input:
        igh_tsv = MIXCR_OUTPUT_DIR + "/{sample}_VH_assemble_Clones_IGH.tsv",
        igk_tsv = MIXCR_OUTPUT_DIR + "/{sample}_VK_assemble_Clones_IGK.tsv",
    output:
        complete_flag = ANA_OUTPUT_DIR + "/{sample}_analysis_complete.txt"
    params:
        data_dir = lambda wildcards: SAMPLES[wildcards.sample],
        output_dir = ANA_OUTPUT_DIR,
        readcount_threshold = config["ana_readcount_threshold"],
        top_n = config["ana_top_n"],
        sample_names = lambda wildcards: [f"{wildcards.sample}-VH", f"{wildcards.sample}-VK"]  # For --samples arg
    log:
        ANA_OUTPUT_DIR + "/logs/{sample}_ana.log"
    shell:
        """
        python ana_modified.py \\
            --data_dir {params.data_dir} \\
            --output_dir {params.output_dir} \\
            --samples {params.sample_names[0]} {params.sample_names[1]} \\
            --readcount_threshold {params.readcount_threshold} \\
            --top_n {params.top_n} &> {log}
        touch {output.complete_flag}
        """

# Rule 3: Aggregate ana summaries (if multi-sample, run once after all)
rule aggregate_ana_summary:
    input:
        expand(ANA_OUTPUT_DIR + "/{sample}_analysis_complete.txt", sample=SAMPLES.keys())
    output:
        summary_xlsx = ANA_OUTPUT_DIR + "/top_combinations_summary.xlsx"
    shell:
        """
        # If ana_modified.py generates per-sample CSVs, merge them here using pandas
        # For simplicity, assume ana_modified.py handles aggregation when --samples includes all
        # Or: python -c "
        # import pandas as pd; from glob import glob;
        # files = glob('{ANA_OUTPUT_DIR}/*.csv'); dfs = [pd.read_csv(f) for f in files];
        # with pd.ExcelWriter('{output.summary_xlsx}') as w: for i,df in enumerate(dfs): df.to_excel(w, sheet_name=f'Sheet{i}', index=False)
        # " 
        echo "Aggregation complete" > {output.summary_xlsx}.flag  # Placeholder; implement merge as needed
        """

# Rule 4: VDJ comparison
rule vdj_compare:
    input:
        unimm_summary = ancient(UNIMMUNIZED_SUMMARY),
        imm_summary = ancient(IMMUNIZED_SUMMARY)
    output:
        top_vh_csv = VDJ_OUTPUT_DIR + "/top_vh_differences.csv",
        top_vk_csv = VDJ_OUTPUT_DIR + "/top_vk_differences.csv"
    params:
        output_dir = VDJ_OUTPUT_DIR,
        top_n = config["vdj_top_n"],
        min_samples = config["vdj_min_samples"],
        alpha = config["vdj_alpha"],
        unimm_vh_sheets = " ".join(UNIMM_VH_SHEETS),
        unimm_vk_sheets = " ".join(UNIMM_VK_SHEETS),
        imm_vh_sheets = " ".join(IMM_VH_SHEETS),
        imm_vk_sheets = " ".join(IMM_VK_SHEETS)
    log:
        VDJ_OUTPUT_DIR + "/logs/vdj_compare.log"
    shell:
        """
        python vdj_modified.py \\
            --unimm_file {input.unimm_summary} \\
            --imm_file {input.imm_summary} \\
            --unimm_vh_sheets {params.unimm_vh_sheets} \\
            --unimm_vk_sheets {params.unimm_vk_sheets} \\
            --imm_vh_sheets {params.imm_vh_sheets} \\
            --imm_vk_sheets {params.imm_vk_sheets} \\
            --top_n {params.top_n} \\
            --min_samples {params.min_samples} \\
            --alpha {params.alpha} \\
            --output_dir {params.output_dir} &> {log}
        """

# Utility rules
rule clean:
    shell:
        """
        rm -rf {MIXCR_OUTPUT_DIR}/logs/ {ANA_OUTPUT_DIR}/logs/ {VDJ_OUTPUT_DIR}/logs/
        find . -name "*.flag" -delete  # Remove temp flags
        """
