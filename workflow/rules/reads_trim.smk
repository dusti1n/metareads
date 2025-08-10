# Rule: [reads_trim_cutadapt]; for illumina and hybrid data_types
rule cutadapt_reads_trim_I_H:
    input:
        r1_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R1.fastq",
        r2_illumina="results_metareads/{project_name}/unpacked/{sample}/{sample}_R2.fastq"
    output:
        trim_r1=temp("results_metareads/{project_name}/filtered/cutadapt/{sample}/{sample}_R1_trim.fastq"),
        trim_r2=temp("results_metareads/{project_name}/filtered/cutadapt/{sample}/{sample}_R2_trim.fastq")
    conda:
        "../envs/reads_trim.yaml"
    threads: config["settings"]["threads"]
    params:
        adapter_r1=config["reads_trim"]["cutadapt"]["adapter_r1"],
        adapter_r2=config["reads_trim"]["cutadapt"]["adapter_r2"],
        min_length=config["reads_trim"]["cutadapt"]["min_length"],
        quality_cutoff=config["reads_trim"]["cutadapt"]["quality_cut_off"],
        overlap=config["reads_trim"]["cutadapt"]["overlap"],
        times=config["reads_trim"]["cutadapt"]["times"]
    shell:
        """
        cutadapt \
            -a {params.adapter_r1} \
            -A {params.adapter_r2} \
            -o {output.trim_r1} \
            -p {output.trim_r2} \
            -q {params.quality_cutoff} \
            --minimum-length {params.min_length} \
            --overlap {params.overlap} \
            --times {params.times} \
            --cores {threads} \
            {input.r1_illumina} \
            {input.r2_illumina} \
        """

# Rule: [porechop_reads_trim_N_H]; for nanopore and hybrid data_types
rule porechop_reads_trim_N_H:
    input:
        ont="results_metareads/{project_name}/unpacked/{sample}/{sample}_{suffix}.fastq"
    output:
        trim_ont=temp("results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_{suffix}_trim.fastq")
    conda:
        "../envs/reads_trim.yaml"
    threads: config["settings"]["threads"] # Threads_(default): porechop=1; no multithread support!
    params:
        adapter_threshold=config["reads_trim"]["porechop"]["adapter_threshold"],
        middle_threshold=config["reads_trim"]["porechop"]["middle_threshold"]
    shell:
        """
        porechop \
            -i {input.ont} \
            -o {output.trim_ont} \
            --threads {threads} \
            --adapter_threshold {params.adapter_threshold} \
            --middle_threshold {params.middle_threshold}
        """
