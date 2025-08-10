# Data type: Nanopore
# Pathway (Flye): Flye-Minimap2-Racon-Medaka

# Rule: [medaka_flye_N]
# Path: [Medaka] after [Racon_Iteration]
rule medaka_flye_N:
    input:
        # Input_file: racon_flye_N_rnd_
        racon_flye_N_rnd_="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_N_rnd_" + str(config["reads_polish"]["racon"]["iterations"]) + ".fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_medaka_N="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_flye_N_final_assemb.fasta"
    conda:
        "../envs/medaka.yaml"
    threads: config["settings"]["threads"] # MEDAKA supports multithreading; threads set in config
    params:
        out_dir="results_metareads/{project_name}/filtered/medaka/{sample}",
    shell:
        """
        export CUDA_VISIBLE_DEVICES=''
        medaka_consensus \
            -i {input.porechop_trim} \
            -d {input.racon_flye_N_rnd_} \
            -o {params.out_dir} \
            -t {threads}
        cp {params.out_dir}/consensus.fasta {output.polished_medaka_N}
        """



# Data type: Nanopore
# Pathway (Canu): Canu-Minimap2-Racon-Medaka

# Rule: [medaka_canu_N]
# Path: [Medaka] after [Racon_Iteration]
rule medaka_canu_N:
    input:
        # Input_file: racon_canu_N_rnd_
        racon_canu_N_rnd_="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_canu_N_rnd_" + str(config["reads_polish"]["racon"]["iterations"]) + ".fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_medaka_N="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_canu_N_final_assemb.fasta"
    conda:
        "../envs/medaka.yaml"
    threads: config["settings"]["threads"] # MEDAKA supports multithreading; threads set in config
    params:
        out_dir="results_metareads/{project_name}/filtered/medaka/{sample}",
    shell:
        """
        export CUDA_VISIBLE_DEVICES=''
        medaka_consensus \
            -i {input.porechop_trim} \
            -d {input.racon_canu_N_rnd_} \
            -o {params.out_dir} \
            -t {threads}
        cp {params.out_dir}/consensus.fasta {output.polished_medaka_N}
        """



# Data type: Hybrid
# Pathway (Flye): Flye-Minimap2-Racon-Medaka

# Rule: [medaka_flye_H]
# Path: [Medaka] after [Racon_Iteration]
rule medaka_flye_H:
    input:
        # Input_file: racon_flye_H_rnd_
        racon_flye_H_rnd_="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_H_rnd_" + str(config["reads_polish"]["racon"]["iterations"]) + ".fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_medaka_H="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_flye_H_assemb.fasta"
    conda:
        "../envs/medaka.yaml"
    threads: config["settings"]["threads"] # MEDAKA supports multithreading; threads set in config
    params:
        out_dir="results_metareads/{project_name}/filtered/medaka/{sample}",
    shell:
        """
        export CUDA_VISIBLE_DEVICES=''
        medaka_consensus \
            -i {input.porechop_trim} \
            -d {input.racon_flye_H_rnd_} \
            -o {params.out_dir} \
            -t {threads}
        cp {params.out_dir}/consensus.fasta {output.polished_medaka_H}
        """
