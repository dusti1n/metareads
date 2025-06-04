# Rule: cutadapt[Illumina, Hybrid]
# Threads: cutadapt=config["settings"]["threads"]
rule cutadapt:
    input:
        r1="results/{project_name}/unpacked/{sample}/{sample}_R1.fastq",
        r2="results/{project_name}/unpacked/{sample}/{sample}_R2.fastq"
    output:
        trimmed_r1="results/{project_name}/filtered/cutadapt/{sample}/{sample}_R1_trimmed.fastq",
        trimmed_r2="results/{project_name}/filtered/cutadapt/{sample}/{sample}_R2_trimmed.fastq"
    conda:
        "../envs/prepolish.yaml"
    threads: config["settings"]["threads"]
    params:
        adapter_r1=config["cutadapt"]["adapter_r1"],
        adapter_r2=config["cutadapt"]["adapter_r2"],
        min_length=config["cutadapt"]["min_length"],
        quality_cutoff=config["cutadapt"]["quality_cut_off"],
        overlap=config["cutadapt"]["overlap"],
        times=config["cutadapt"]["times"]
    shell:
        """
        cutadapt \
            -a {params.adapter_r1} \
            -A {params.adapter_r2} \
            -o {output.trimmed_r1} \
            -p {output.trimmed_r2} \
            -q {params.quality_cutoff} \
            --minimum-length {params.min_length} \
            --overlap {params.overlap} \
            --times {params.times} \
            --cores {threads} \
            {input.r1} \
            {input.r2} \
        """

# Rule: porechop[Nanopore, Hybrid]
# Threads (default): porechop=1
rule porechop:
    input:
        ont="results/{project_name}/unpacked/{sample}/{sample}_{suffix}.fastq"
    output:
        trimmed_ont="results/{project_name}/filtered/porechop/{sample}/{sample}_{suffix}_trimmed.fastq"
    conda:
        "../envs/prepolish.yaml"
    threads: config["settings"]["threads"]
    params:
        adapter_threshold=config["porechop"]["adapter_threshold"],
        middle_threshold=config["porechop"]["middle_threshold"]
    shell:
        """
        porechop \
            -i {input.ont} \
            -o {output.trimmed_ont} \
            --threads {threads} \
            --adapter_threshold {params.adapter_threshold} \
            --middle_threshold {params.middle_threshold}
        """
