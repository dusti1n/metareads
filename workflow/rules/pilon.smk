# Data type: Illumina
# Pathway (Spades): Spades-BWA-Pilon

# Rule: [bwa_map_spades_I]
# Path: [bwa_map_spades_I] before [pilon_spades_I]
rule bwa_map_spades_I:
    input:
        # Input_file: assembly_spades
        assembly_spades="results_metareads/{project_name}/assembled/spades/{sample}/{sample}_assemb.fasta",
        # Input_file: r1_illumina, r2_illumina
        r1_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R1.fastq",
        r2_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R2.fastq"
    output:
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_spades_I.bam"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # BWA supports multithreading; threads set in config
    shell:
        """
        bwa index {input.assembly_spades}
        bwa mem -M -t {threads} {input.assembly_spades} {input.r1_illumina} {input.r2_illumina} | \
            samtools view -Sb -u - | \
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# Rule: [pilon_spades_I]
# Path: [pilon_spades_I] after [bwa_map_spades_I] 
rule pilon_spades_I:
    input:
        # Input_file: bam file from rule bwa_map_spades_I
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_spades_I.bam",
        # Input_file: assembly_spades
        assembly_spades="results_metareads/{project_name}/assembled/spades/{sample}/{sample}_assemb.fasta"
    output:
        polished_spades_I="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_pilon_spades_I_final_assemb.fasta"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # Still used for Snakemake resource management
    params:
        outdir="results_metareads/{project_name}/filtered/pilon/{sample}",
        prefix="{sample}_pilon_spades_I_final_assemb",
        heap=int(config["settings"]["memory"]) // 1000,
        fix="--fix " + config["reads_polish"]["pilon"]["fix"] if config["reads_polish"]["pilon"].get("fix") else "",
        diploid="--diploid" if config["reads_polish"]["pilon"].get("diploid") else "",
        mindepth="--mindepth " + str(config["reads_polish"]["pilon"]["mindepth"]) if config["reads_polish"]["pilon"].get("mindepth") is not None else "",
        minmq="--minmq " + str(config["reads_polish"]["pilon"]["minmq"]) if config["reads_polish"]["pilon"].get("minmq") is not None else "",
        minqual="--minqual " + str(config["reads_polish"]["pilon"]["minqual"]) if config["reads_polish"]["pilon"].get("minqual") is not None else "",
        chunksize="--chunksize " + str(config["reads_polish"]["pilon"]["chunksize"]) if config["reads_polish"]["pilon"].get("chunksize") is not None else ""
    shell:
        """
        # Index BAM file
        samtools index {input.bam}

        # Set Java heap size
        export _JAVA_OPTIONS="-Xmx{params.heap}G"

        # Run Pilon for assembly polishing
        pilon \
        --genome {input.assembly_spades} \
        --bam {input.bam} \
        --outdir {params.outdir} \
        --output {params.prefix} \
        {params.fix} \
        {params.diploid} \
        {params.mindepth} \
        {params.minmq} \
        {params.minqual} \
        {params.chunksize} \
        --changes \
        --vcf
        """



# Data type: Illumina
# Pathway (MEGAHIT): MEGAHIT-BWA-Pilon

# Rule: [bwa_map_megahit_I]
# Path: [bwa_map_megahit_I] before [pilon_megahit_I]
rule bwa_map_megahit_I:
    input:
        # Input_file: assembly_megahit
        assembly_megahit="results_metareads/{project_name}/assembled/megahit/{sample}/{sample}_assemb.contigs.fa",
        # Input_file: r1_illumina, r2_illumina
        r1_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R1.fastq",
        r2_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R2.fastq"
    output:
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_megahit_I.bam"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # BWA supports multithreading; threads set in config
    shell:
        """
        bwa index {input.assembly_megahit}
        bwa mem -M -t {threads} {input.assembly_megahit} {input.r1_illumina} {input.r2_illumina} | \
            samtools view -Sb -u - | \
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# Rule: [pilon_megahit_I]
# Path: [pilon_megahit_I] after [bwa_map_megahit_I] 
rule pilon_megahit_I:
    input:
        # Input_file: bam file from rule bwa_map_megahit_I
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_megahit_I.bam",
        # Input_file: assembly_megahit
        assembly_megahit="results_metareads/{project_name}/assembled/megahit/{sample}/{sample}_assemb.contigs.fa"
    output:
        polished_megahit_I="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_pilon_megahit_I_final_assemb.fasta"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # Still used for Snakemake resource management
    params:
        outdir="results_metareads/{project_name}/filtered/pilon/{sample}",
        prefix="{sample}_pilon_megahit_I_final_assemb",
        heap=int(config["settings"]["memory"]) // 1000,
        fix="--fix " + config["reads_polish"]["pilon"]["fix"] if config["reads_polish"]["pilon"].get("fix") else "",
        diploid="--diploid" if config["reads_polish"]["pilon"].get("diploid") else "",
        mindepth="--mindepth " + str(config["reads_polish"]["pilon"]["mindepth"]) if config["reads_polish"]["pilon"].get("mindepth") is not None else "",
        minmq="--minmq " + str(config["reads_polish"]["pilon"]["minmq"]) if config["reads_polish"]["pilon"].get("minmq") is not None else "",
        minqual="--minqual " + str(config["reads_polish"]["pilon"]["minqual"]) if config["reads_polish"]["pilon"].get("minqual") is not None else "",
        chunksize="--chunksize " + str(config["reads_polish"]["pilon"]["chunksize"]) if config["reads_polish"]["pilon"].get("chunksize") is not None else ""
    shell:
        """
        # Index the BAM file to enable fast random access
        samtools index {input.bam}

        # Set Java heap size
        export _JAVA_OPTIONS="-Xmx{params.heap}G"

        # Run Pilon for assembly polishing
        pilon \
        --genome {input.assembly_megahit} \
        --bam {input.bam} \
        --outdir {params.outdir} \
        --output {params.prefix} \
        {params.fix} \
        {params.diploid} \
        {params.mindepth} \
        {params.minmq} \
        {params.minqual} \
        {params.chunksize} \
        --changes \
        --vcf
        """



# Data type: Hybrid
# Pathway (hSpades): hSpades-BWA-Pilon

# Rule: [bwa_map_hSpades_H]
# Path: [bwa_map_hSpades_H] before [pilon_hSpades_H]
rule bwa_map_hSpades_H:
    input:
        # Input_file: assembly_hSpades
        assembly_hSpades="results_metareads/{project_name}/assembled/spades/{sample}/{sample}_assemb.fasta",
        # Input_file: r1_illumina, r2_illumina
        r1_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R1.fastq",
        r2_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R2.fastq"
    output:
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_hSpades_H.bam"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # BWA supports multithreading; threads set in config
    shell:
        """
        bwa index {input.assembly_hSpades}
        bwa mem -M -t {threads} {input.assembly_hSpades} {input.r1_illumina} {input.r2_illumina} | \
            samtools view -Sb -u - | \
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# Rule: [pilon_hSpades_H]
# Path: [pilon_hSpades_H] after [bwa_map_hSpades_H] 
rule pilon_hSpades_H:
    input:
        # Input_file: bam file from rule bwa_map_hSpades_H
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_hSpades_H.bam",
        # Input_file: assembly_hSpades
        assembly_hSpades="results_metareads/{project_name}/assembled/spades/{sample}/{sample}_assemb.fasta"
    output:
        polished_hSpades_H="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_hSpades_H_final_assemb.fasta"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # Still used for Snakemake resource management
    params:
        outdir="results_metareads/{project_name}/filtered/pilon/{sample}",
        prefix="{sample}_hSpades_H_final_assemb",
        heap=int(config["settings"]["memory"]) // 1000,
        fix="--fix " + config["reads_polish"]["pilon"]["fix"] if config["reads_polish"]["pilon"].get("fix") else "",
        diploid="--diploid" if config["reads_polish"]["pilon"].get("diploid") else "",
        mindepth="--mindepth " + str(config["reads_polish"]["pilon"]["mindepth"]) if config["reads_polish"]["pilon"].get("mindepth") is not None else "",
        minmq="--minmq " + str(config["reads_polish"]["pilon"]["minmq"]) if config["reads_polish"]["pilon"].get("minmq") is not None else "",
        minqual="--minqual " + str(config["reads_polish"]["pilon"]["minqual"]) if config["reads_polish"]["pilon"].get("minqual") is not None else "",
        chunksize="--chunksize " + str(config["reads_polish"]["pilon"]["chunksize"]) if config["reads_polish"]["pilon"].get("chunksize") is not None else ""
    shell:
        """
        # Index the BAM file to enable fast random access
        samtools index {input.bam}

        # Set Java heap size
        export _JAVA_OPTIONS="-Xmx{params.heap}G"

        # Run Pilon for assembly polishing
        pilon \
        --genome {input.assembly_hSpades} \
        --bam {input.bam} \
        --outdir {params.outdir} \
        --output {params.prefix} \
        {params.fix} \
        {params.diploid} \
        {params.mindepth} \
        {params.minmq} \
        {params.minqual} \
        {params.chunksize} \
        --changes \
        --vcf
        """



# Data type: Hybrid
# Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-(BWA-Pilon)

# Rule: [bwa_map_flye_H]
# Path: [bwa_map_flye_H] before [pilon_flye_H]
rule bwa_map_flye_H:
    input:
        # Input_file: assembly_flye_medaka_H
        assembly_flye_medaka_H="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_flye_H_assemb.fasta",
        # Input_file: r1_illumina, r2_illumina
        r1_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R1.fastq",
        r2_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R2.fastq"
    output:
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_flye_H.bam"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # BWA supports multithreading; threads set in config
    shell:
        """
        bwa index {input.assembly_flye_medaka_H}
        bwa mem -M -t {threads} {input.assembly_flye_medaka_H} {input.r1_illumina} {input.r2_illumina} | \
            samtools view -Sb -u - | \
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# Rule: [pilon_flye_H]
# Path: [pilon_flye_H] after [bwa_map_flye_H] 
rule pilon_flye_H:
    input:
        # Input_file: bam file from rule bwa_map_flye_H
        bam="results_metareads/{project_name}/filtered/bwa_map/{sample}/{sample}_bwa_map_flye_H.bam",
        # Input_file: assembly_flye_medaka_H
        assembly_flye_medaka_H="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_flye_H_assemb.fasta",
    output:
        polished_flye_H="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_flye_H_final_assemb.fasta"
    conda:
        "../envs/pilon.yaml"
    threads: config["settings"]["threads"] # Still used for Snakemake resource management
    params:
        outdir="results_metareads/{project_name}/filtered/pilon/{sample}",
        prefix="{sample}_flye_H_final_assemb",
        heap=int(config["settings"]["memory"]) // 1000,
        fix="--fix " + config["reads_polish"]["pilon"]["fix"] if config["reads_polish"]["pilon"].get("fix") else "",
        diploid="--diploid" if config["reads_polish"]["pilon"].get("diploid") else "",
        mindepth="--mindepth " + str(config["reads_polish"]["pilon"]["mindepth"]) if config["reads_polish"]["pilon"].get("mindepth") is not None else "",
        minmq="--minmq " + str(config["reads_polish"]["pilon"]["minmq"]) if config["reads_polish"]["pilon"].get("minmq") is not None else "",
        minqual="--minqual " + str(config["reads_polish"]["pilon"]["minqual"]) if config["reads_polish"]["pilon"].get("minqual") is not None else "",
        chunksize="--chunksize " + str(config["reads_polish"]["pilon"]["chunksize"]) if config["reads_polish"]["pilon"].get("chunksize") is not None else ""
    shell:
        """
        # Index the BAM file to enable fast random access
        samtools index {input.bam}

        # Set Java heap size
        export _JAVA_OPTIONS="-Xmx{params.heap}G"

        # Run Pilon for assembly polishing
        pilon \
        --genome {input.assembly_flye_medaka_H} \
        --bam {input.bam} \
        --outdir {params.outdir} \
        --output {params.prefix} \
        {params.fix} \
        {params.diploid} \
        {params.mindepth} \
        {params.minmq} \
        {params.minqual} \
        {params.chunksize} \
        --changes \
        --vcf
        """
