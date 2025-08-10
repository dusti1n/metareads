# Rule: [unzip_samples]; for all data_types
rule unzip_samples:
    priority: 500
    input:
        # Dependency: enforces DAG plot generation before this rule
        # dag_plot=os.path.join("results_metareads", config["general"]["output_dir"], "dag_plot_tree.pdf"),
        # Input: compressed FASTQ file
        zipped=lambda wildcards: sample_files[wildcards.sample].get(
            {"R1": "fq1", "R2": "fq2", "ONT": "ONT"}[wildcards.suffix])
    output:
        unzipped=temp("results_metareads/{project_name}/unpacked/{sample}/{sample}_{suffix}.fastq")
    conda:
        "../envs/unzip_samples.yaml"
    threads: config["settings"]["threads"]
    shell:
        """
        pigz -p {threads} -d -c {input.zipped} > {output.unzipped}
        """
