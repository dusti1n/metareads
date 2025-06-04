# Rule: asmbflye[Nanopore, Hybrid]
rule asmbflye:
    input:
        trimmed_ont="results/{project_name}/filtered/porechop/{sample}/{sample}_{suffix}_trimmed.fastq"
    output:
        assembly_flye="results/{project_name}/assembles/flye/{sample}/{sample}_{suffix}_asmb.fasta"
    conda:
        "../envs/asmbflye.yaml"
    params:
        flye_dir="results/{project_name}/assembles/flye/{sample}/",
        genome_size=config["assembler"]["flye"]["genome_size"],
        iterations=config["assembler"]["flye"]["iterations"]
    threads: config["settings"]["threads"]
    shell:
        """
        flye --nano-raw {input.trimmed_ont} \
             --out-dir {params.flye_dir} \
             --genome-size {params.genome_size} \
             --iterations {params.iterations} \
             --threads {threads}
        
        mv {params.flye_dir}/assembly.fasta {output.assembly_flye}
        """

# Rule: asmbcanu[Nanopore, Hybrid]
rule asmbcanu:
    input:
        trimmed_ont="results/{project_name}/filtered/porechop/{sample}/{sample}_{suffix}_trimmed.fastq"
    output:
        assembly_canu="results/{project_name}/assembles/canu/{sample}/{sample}_{suffix}.contigs.fasta"
    conda:
        "../envs/asmbcanu.yaml"
    params:
        canu_dir="results/{project_name}/assembles/canu/{sample}/",
        prefix="{sample}_{suffix}",
        genome_size=config["assembler"]["canu"]["genome_size"],
        low_coverage=config["assembler"]["canu"]["stop_on_low_coverage"],
        min_coverage=config["assembler"]["canu"]["min_input_coverage"],
    threads: config["settings"]["threads"]
    shell:
        """
        canu -p {params.prefix} -d {params.canu_dir} \
            genomeSize={params.genome_size} \
            -nanopore {input.trimmed_ont} \
            useGrid=false \
            stopOnLowCoverage={params.low_coverage} \
            minInputCoverage={params.min_coverage}
        """

# Rule: asmbspades[Illumina, Hybrid]
rule asmbspades:
    input:
        trimmed_r1="results/{project_name}/filtered/cutadapt/{sample}/{sample}_R1_trimmed.fastq",
        trimmed_r2="results/{project_name}/filtered/cutadapt/{sample}/{sample}_R2_trimmed.fastq",
        trimmed_ont="results/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trimmed.fastq" if config["settings"]["data_type"] == "Hybrid" else[]
    output:
        assembly_spades="results/{project_name}/assembles/spades/{sample}/{sample}_asmb.fasta"
    conda:
        "../envs/asmbspades.yaml"
    params:
        data_type=config["settings"]["data_type"],
        spades_dir="results/{project_name}/assembles/spades/{sample}",
    threads: config["settings"]["threads"]
    shell:
        """
        mkdir -p {params.spades_dir}

        if [[ "{params.data_type}" == "Hybrid" && -s "{input.trimmed_ont}" ]]; then
            echo "INFO: Running HybridSPAdes (Illumina + Nanopore)"
            spades.py --nanopore {input.trimmed_ont} \
                      -1 {input.trimmed_r1} -2 {input.trimmed_r2} \
                      --isolate -o {params.spades_dir} --threads {threads}
        else
            echo "INFO: Running SPAdes in Illumina-Only Mode"
            spades.py -1 {input.trimmed_r1} -2 {input.trimmed_r2} \
                      --isolate -o {params.spades_dir} --threads {threads}
        fi

        mv {params.spades_dir}/contigs.fasta {output.assembly_spades}
        """

# Rule: asmbmegahit[Illumina]
rule asmbmegahit:
    input:
        trimmed_r1="results/{project_name}/filtered/cutadapt/{sample}/{sample}_R1_trimmed.fastq",
        trimmed_r2="results/{project_name}/filtered/cutadapt/{sample}/{sample}_R2_trimmed.fastq",
    output:
        assembly_megahit="results/{project_name}/assembles/megahit/{sample}/{sample}_asmb.contigs.fa"
    conda:
        "../envs/asmbmegahit.yaml"
    params:
        out_dir="results/{project_name}/assembles/megahit/{sample}",
        memory=config["settings"]["memory"],
        min_contig_len=config["assembler"]["megahit"]["min_contig_len"],
        min_count=config["assembler"]["megahit"]["min_count"],
        k_min=config["assembler"]["megahit"]["k_min"],
        k_max=config["assembler"]["megahit"]["k_max"],
        k_step=config["assembler"]["megahit"]["k_step"]
    threads: config["settings"]["threads"]
    shell:
        """
        rm -rf {params.out_dir}

        megahit -1 {input.trimmed_r1} \
                -2 {input.trimmed_r2} \
                --out-dir {params.out_dir} \
                --num-cpu-threads {threads} \
                --memory {params.memory} \
                --min-contig-len {params.min_contig_len} \
                --min-count {params.min_count} \
                --k-min {params.k_min} \
                --k-max {params.k_max} \
                --k-step {params.k_step}

        mv {params.out_dir}/final.contigs.fa {output.assembly_megahit}
        """
