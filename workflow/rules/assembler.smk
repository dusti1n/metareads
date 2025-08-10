# Rule: [assembler_flye_N_H]; N=Nanopore; H=Hybrid; for nanopore and hybrid data_types
rule assembler_flye_N_H:
    input:
        trim_ont="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_{suffix}_trim.fastq"
    output:
        assembly_flye="results_metareads/{project_name}/assembled/flye/{sample}/{sample}_{suffix}_assemb.fasta"
    conda:
        "../envs/assembler_flye.yaml"
    params:
        flye_dir="results_metareads/{project_name}/assembled/flye/{sample}/",
        genome_size=config["assembler"]["flye"]["genome_size"],
        iterations=config["assembler"]["flye"]["iterations"],
        meta_flag="--meta" if config["assembler"]["flye"].get("meta", False) else ""
    threads: config["settings"]["threads"]
    shell:
        """
        flye --nano-raw {input.trim_ont} \
             --out-dir {params.flye_dir} \
             --genome-size {params.genome_size} \
             --iterations {params.iterations} \
             --threads {threads} \
             {params.meta_flag}
        
        mv {params.flye_dir}/assembly.fasta {output.assembly_flye}
        """

# Rule: [assembler_canu_N]; N=Nanopore; for nanopore data_type
rule assembler_canu_N:
    input:
        trim_ont="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_{suffix}_trim.fastq"
    output:
        assembly_canu="results_metareads/{project_name}/assembled/canu/{sample}/{sample}_{suffix}_assemb.contigs.fasta"
    conda:
        "../envs/assembler_canu.yaml"
    params:
        canu_dir="results_metareads/{project_name}/assembled/canu/{sample}/",
        prefix="{sample}_{suffix}_assemb",
        genome_size=config["assembler"]["canu"]["genome_size"],
        low_coverage=config["assembler"]["canu"]["stop_on_low_coverage"],
        min_coverage=config["assembler"]["canu"]["min_input_coverage"],
    threads: config["settings"]["threads"]
    shell:
        """
        canu -p {params.prefix} -d {params.canu_dir} \
            genomeSize={params.genome_size} \
            -nanopore {input.trim_ont} \
            useGrid=false \
            stopOnLowCoverage={params.low_coverage} \
            minInputCoverage={params.min_coverage}
        """

# Rule: [assembler_spades_I_H]; I=Illumina; H=Hybrid; for illumina and hybrid data_types
rule assembler_spades_I_H:
    input:
        trim_r1="results_metareads/{project_name}/filtered/cutadapt/{sample}/{sample}_R1_trim.fastq",
        trim_r2="results_metareads/{project_name}/filtered/cutadapt/{sample}/{sample}_R2_trim.fastq",
        trim_ont="results_metareads/{project_name}/filtered/porechop/{sample}/{sample}_ONT_trim.fastq" if config["settings"]["data_type"] == "hybrid" else[]
    output:
        assembly_spades="results_metareads/{project_name}/assembled/spades/{sample}/{sample}_assemb.fasta"
    conda:
        "../envs/assembler_spades.yaml"
    params:
        data_type=config["settings"]["data_type"],
        spades_dir="results_metareads/{project_name}/assembled/spades/{sample}",
        meta_mode="--meta" if config["assembler"]["spades"]["meta"] else "--isolate"
    threads: config["settings"]["threads"]
    shell:
        """
        mkdir -p {params.spades_dir}

        if [[ "{params.data_type}" == "hybrid" && -s "{input.trim_ont}" ]]; then
            echo "INFO: Running HybridSPAdes (Illumina with Nanopore)"
            spades.py --nanopore {input.trim_ont} \
                      -1 {input.trim_r1} -2 {input.trim_r2} \
                      -o {params.spades_dir} --threads {threads}
        else
            echo "INFO: Running SPAdes in Illumina-Only Mode"
            spades.py -1 {input.trim_r1} -2 {input.trim_r2} \
                      {params.meta_mode} \
                      -o {params.spades_dir} --threads {threads}
        fi

        mv {params.spades_dir}/contigs.fasta {output.assembly_spades}
        """

# Rule: [assembler_megahit_I]; I=Illumina; for illumina data_type
rule assembler_megahit_I:
    input:
        trim_r1="results_metareads/{project_name}/filtered/cutadapt/{sample}/{sample}_R1_trim.fastq",
        trim_r2="results_metareads/{project_name}/filtered/cutadapt/{sample}/{sample}_R2_trim.fastq"
    output:
        assembly_megahit="results_metareads/{project_name}/assembled/megahit/{sample}/{sample}_assemb.contigs.fa"
    conda:
        "../envs/assembler_megahit.yaml"
    params:
        out_dir="results_metareads/{project_name}/assembled/megahit/{sample}",
        memory=config["settings"]["memory"],
        min_contig_len="--min-contig-len " + str(config["assembler"]["megahit"]["min_contig_len"]),
        megahit_flags=(
            "--min-count {min_count} --k-min {k_min} --k-max {k_max} --k-step {k_step}".format(
                min_count=config["assembler"]["megahit"]["min_count"],
                k_min=config["assembler"]["megahit"]["k_min"],
                k_max=config["assembler"]["megahit"]["k_max"],
                k_step=config["assembler"]["megahit"]["k_step"]
            )
            if config["assembler"]["megahit"]["preset"] == "manual"
            else "--preset " + config["assembler"]["megahit"]["preset"]
        )
    threads: config["settings"]["threads"]
    shell:
        """
        rm -rf {params.out_dir}

        megahit -1 {input.trim_r1} \
                -2 {input.trim_r2} \
                --out-dir {params.out_dir} \
                --num-cpu-threads {threads} \
                --memory {params.memory} \
                {params.min_contig_len} \
                {params.megahit_flags}

        mv {params.out_dir}/final.contigs.fa {output.assembly_megahit}
        """
