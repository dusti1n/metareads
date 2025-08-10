# Data type: Nanopore
# Pathway (Flye): Flye-Minimap2-Racon
# Iteration: Racon_1-3

# Rule: [racon_flye_N_rnd_1]
# Iteration_1; [racon_flye_N_rnd_1] after [flye_assembler]
rule racon_flye_N_rnd_1:
    input:
        # Input_file: flye_assembly
        flye_assembly="results_metareads/{project_name}/assembled/flye/{sample}/{sample}_ONT_assemb.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_N_rnd_1.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_flye_N_rnd_1_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.flye_assembly} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.flye_assembly} > {output.polished_racon}
        """

# Rule: [racon_flye_N_rnd_2] 
# Iteration_2; [racon_flye_N_rnd_2] after [racon_flye_N_rnd_1]
rule racon_flye_N_rnd_2:
    input:
        # Input_file: racon_flye_N_rnd_1
        racon_flye_N_rnd_1="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_N_rnd_1.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_N_rnd_2.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_flye_N_rnd_2_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.racon_flye_N_rnd_1} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.racon_flye_N_rnd_1} > {output.polished_racon}
        """

# Rule: [racon_flye_N_rnd_3]
# Iteration_3; [racon_flye_N_rnd_3] after [racon_flye_N_rnd_2]
rule racon_flye_N_rnd_3:
    input:
        # Input_file: racon_flye_N_rnd_2
        racon_flye_N_rnd_2="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_N_rnd_2.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_N_rnd_3.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_flye_N_rnd_3_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.racon_flye_N_rnd_2} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.racon_flye_N_rnd_2} > {output.polished_racon}
        """



# Data type: Nanopore
# Pathway (Canu): Canu-Minimap2-Racon
# Iteration: Racon_1-3

# Rule: [racon_canu_N_rnd_1]
# Iteration_1; [racon_canu_N_rnd_1] after [canu_assembly]
rule racon_canu_N_rnd_1:
    input:
        # Input_file: canu_assembly
        canu_assembly="results_metareads/{project_name}/assembled/canu/{sample}/{sample}_ONT_assemb.contigs.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_canu_N_rnd_1.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_canu_N_rnd_1_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.canu_assembly} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.canu_assembly} > {output.polished_racon}
        """

# Rule: [racon_canu_N_rnd_2] 
# Iteration_2; [racon_canu_N_rnd_2] after [racon_canu_N_rnd_1]
rule racon_canu_N_rnd_2:
    input:
        # Input_file: racon_canu_N_rnd_1
        racon_canu_N_rnd_1="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_canu_N_rnd_1.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_canu_N_rnd_2.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_canu_N_rnd_2_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.racon_canu_N_rnd_1} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.racon_canu_N_rnd_1} > {output.polished_racon}
        """

# Rule: [racon_canu_N_rnd_3]
# Iteration_3; [racon_canu_N_rnd_3] after [racon_canu_N_rnd_2]
rule racon_canu_N_rnd_3:
    input:
        # Input_file: racon_canu_N_rnd_2
        racon_canu_N_rnd_2="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_canu_N_rnd_2.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_canu_N_rnd_3.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_canu_N_rnd_3_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.racon_canu_N_rnd_2} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.racon_canu_N_rnd_2} > {output.polished_racon}
        """



# Data type: Hybrid
# Pathway (Flye): Flye-Minimap2-Racon
# Iteration: Racon_1-3

# Rule: [racon_flye_H_rnd_1]
# Iteration_1; [racon_flye_H_rnd_1] after [flye]
rule racon_flye_H_rnd_1:
    input:
        # Input_file: flye_assembly
        flye_assembly="results_metareads/{project_name}/assembled/flye/{sample}/{sample}_ONT_assemb.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_H_rnd_1.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_flye_H_rnd_1_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.flye_assembly} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.flye_assembly} > {output.polished_racon}
        """

# Rule: [racon_flye_H_rnd_2]
# Iteration_2; [racon_flye_H_rnd_2] after [racon_flye_H_rnd_1]
rule racon_flye_H_rnd_2:
    input:
        # Input_file: racon_flye_H_rnd_1
        racon_flye_H_rnd_1="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_H_rnd_1.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_H_rnd_2.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_flye_H_rnd_2_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.racon_flye_H_rnd_1} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.racon_flye_H_rnd_1} > {output.polished_racon}
        """

# Rule: [racon_flye_H_rnd_3]
# Iteration_3; [racon_flye_H_rnd_3] after [racon_flye_H_rnd_2]
rule racon_flye_H_rnd_3:
    input:
        # Input_file: racon_flye_H_rnd_2
        racon_flye_H_rnd_2="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_H_rnd_2.fasta",
        # Input_file: porechop_trim
        porechop_trim="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq"
    output:
        polished_racon="results_metareads/{project_name}/filtered/racon/{sample}/{sample}_racon_flye_H_rnd_3.fasta",
        sam=temp("results_metareads/{project_name}/filtered/racon/{sample}_racon_flye_H_rnd_3_map.sam"), # Used for Racon polishing (input alignment)
    conda:
        "../envs/racon.yaml"
    threads: config["settings"]["threads"] # Racon supports multithreading; threads set in config
    params:
        match=config["reads_polish"]["racon"]["match_score"],
        mismatch=config["reads_polish"]["racon"]["mismatch_penalty"],
        quality = config["reads_polish"]["racon"]["quality_threshold"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.racon_flye_H_rnd_2} {input.porechop_trim} > {output.sam}
        racon -t {threads} \
              -m {params.match} \
              -x {params.mismatch} \
              -q {params.quality} \
              {input.porechop_trim} {output.sam} {input.racon_flye_H_rnd_2} > {output.polished_racon}
        """
