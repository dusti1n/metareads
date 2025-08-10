# Rule: [fastqc_reads_quality_I_H]; for illumina and hybrid data_types
# Threads_(default): fastqc=1
rule fastqc_reads_quality_I_H:
    input:
        fastq="results_metareads/{project_name}/unpacked/{sample}/{sample}_{suffix}.fastq"
    output:
        html=report("results_metareads/{project_name}/quality/{sample}/fastqc/{sample}_{suffix}_fastqc.html", category="FastQC")
    conda:
        "../envs/reads_quality.yaml"
    params:
        outdir="results_metareads/{project_name}/quality/{sample}/fastqc/"
    shell:
        """
        fastqc {input.fastq} --outdir={params.outdir}
        """

# Rule: [multiqc_reads_quality_I_H]; for illumina and hybrid data_types
# Threads_(default): multiqc=1
rule multiqc_reads_quality_I_H:
    input:
        expand("results_metareads/{project_name}/quality/{sample}/fastqc/{sample}_{suffix}_fastqc.html",
               sample=sample_files.keys(),
               project_name=config["general"]["output_dir"],
               suffix=["R1", "R2"]) if config["settings"]['data_type'] in ["illumina", "hybrid"] else[]
    output:
        report(os.path.join("results_metareads",config["general"]["output_dir"],"quality/multiqc/multiqc_report.html"), category="MultiQC")
    conda:
        "../envs/reads_quality.yaml"
    params:
        project_name=config["general"]["output_dir"]
    shell:
        """
        multiqc results_metareads/{params.project_name}/quality/ -o results_metareads/{params.project_name}/quality/multiqc/
        """

# Rule: [nanoplot_reads_quality_N_H]; for nanopore and hybrid data_types
# Threads_(default): nanoplot=4
rule nanoplot_reads_quality_N_H:
    input:
        fastq="results_metareads/{project_name}/unpacked/{sample}/{sample}_{suffix}.fastq"
    output:
        report=report("results_metareads/{project_name}/quality/{sample}/nanoplot/{sample}_{suffix}_nanoplot.html", category="Nanoplot")
    conda:
        "../envs/reads_quality.yaml"
    params:
        outdir="results_metareads/{project_name}/quality/{sample}/nanoplot/",
        prefix="{sample}_{suffix}_nanoplot"
    shell:
        """
        NanoPlot --fastq {input.fastq} \
        --outdir={params.outdir} \
        --prefix {params.prefix}
        mv {params.outdir}/{params.prefix}NanoPlot-report.html {output.report}
        """
