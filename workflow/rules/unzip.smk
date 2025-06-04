# Rule: unzip[Illumina, Nanopore, Hybrid]
# Threads: pigz=config["settings"]["threads"]
rule unzip:
    input:
        gzipped=lambda wildcards: sample_files[wildcards.sample].get(
            {"R1": "fq1", "R2": "fq2", "ONT": "ONT"}[wildcards.suffix]
        )
    output:
        unzipped="results/{project_name}/unpacked/{sample}/{sample}_{suffix}.fastq"
    conda:
        "../envs/unzip.yaml"
    threads: config["settings"]["threads"]
    shell:
        """
        pigz -p {threads} -d -c {input.gzipped} > {output.unzipped}
        """
