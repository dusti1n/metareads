# Data type: Illumina
# Pathway (Spades/MEGAHIT): Spades/MEGAHIT-BWA-Pilon-(QUAST)
# Analyze final assemblies with QUAST

if config["settings"]["data_type"] == "illumina":
    # Rule: [quast_analysis_I]; for illumina data type
    rule quast_analysis_I:
        input:
            # Pathway (Spades): Spades-BWA-Pilon-(QUAST)
            # Input_file: {sample}_pilon_spades_I_final_assemb.fasta
            polished_spades_I="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_pilon_spades_I_final_assemb.fasta"
                if config["settings"]["data_type"] == "illumina" and config["assembler"]["spades"]["status"] else [],

            # Pathway (MEGAHIT): MEGAHIT-BWA-Pilon-(QUAST)
            # Input_file: {sample}_pilon_megahit_I_final_assemb.fasta
            polished_megahit_I="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_pilon_megahit_I_final_assemb.fasta"
                if config["settings"]["data_type"] == "illumina" and config["assembler"]["megahit"]["status"] else []
        output:
            html=report("results_metareads/{project_name}/analysis/quast/{sample}_quast_eval.html", category="QUAST"),
            tsv=report("results_metareads/{project_name}/analysis/quast/{sample}_quast_eval.tsv", category="QUAST"),
        conda:
            "../envs/quast.yaml"
        threads: config["settings"]["threads"]
        params:
            outdir="results_metareads/{project_name}/analysis/quast"
        shell:
            """
            quast {input.polished_spades_I} {input.polished_megahit_I} -o {params.outdir} -t {threads}
            mv {params.outdir}/report.html {output.html}
            mv {params.outdir}/report.tsv {output.tsv}
            """


# Data type: Nanopore
# Pathway (Flye/Canu): Flye/Canu-Minimap2-Racon-MEDAKA-(QUAST)
# Analyze final assemblies with QUAST

if config["settings"]["data_type"] == "nanopore":
    # Rule: [quast_analysis_N]; for nanopore data type
    rule quast_analysis_N:
        input:
            # Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-(QUAST)
            # Input_file: {sample}_flye_N_final_assemb.fasta
            polished_flye_medaka_N="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_flye_N_final_assemb.fasta"
                if config["settings"]["data_type"] == "nanopore" and config["assembler"]["flye"]["status"] else [],

            # Pathway (Canu): Canu-Minimap2-Racon-MEDAKA-(QUAST)
            # Input_file: {sample}_canu_N_final_assemb.fasta
            polished_canu_medaka_N="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_canu_N_final_assemb.fasta"
                if config["settings"]["data_type"] == "nanopore" and config["assembler"]["canu"]["status"] else []
        output:
            html=report("results_metareads/{project_name}/analysis/quast/{sample}_quast_eval.html", category="QUAST"),
            tsv=report("results_metareads/{project_name}/analysis/quast/{sample}_quast_eval.tsv", category="QUAST"),
        conda:
            "../envs/quast.yaml"
        threads: config["settings"]["threads"]
        params:
            outdir="results_metareads/{project_name}/analysis/quast"
        shell:
            """
            quast {input.polished_flye_medaka_N} {input.polished_canu_medaka_N} -o {params.outdir} -t {threads}
            mv {params.outdir}/report.html {output.html}
            mv {params.outdir}/report.tsv {output.tsv}
            """


# Data type: Hybrid (assemblers: Flye/hSpades)
# Analyze final assemblies with QUAST

if config["settings"]["data_type"] == "hybrid":
    # Rule: [quast_analysis_H]; for hybrid data type
    rule quast_analysis_H:
        input:
            # Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-BWA-Pilon-(QUAST)
            # Input_file: {sample}_flye_H_final_assemb.fasta
            polished_flye_H="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_flye_H_final_assemb.fasta"
                if config["settings"]["data_type"] == "hybrid" and config["assembler"]["flye"]["status"] else [],

            # Pathway (hSpades): hSpades-BWA-Pilon-(QUAST)
            # Input_file: {sample}_hSpades_H_final_assemb.fasta
            polished_hSpades_H="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_hSpades_H_final_assemb.fasta"
                if config["settings"]["data_type"] == "hybrid" and config["assembler"]["spades"]["status"] else []
        output:
            html=report("results_metareads/{project_name}/analysis/quast/{sample}_quast_eval.html", category="QUAST"),
            tsv=report("results_metareads/{project_name}/analysis/quast/{sample}_quast_eval.tsv", category="QUAST"),
        conda:
            "../envs/quast.yaml"
        threads: config["settings"]["threads"]
        params:
            outdir="results_metareads/{project_name}/analysis/quast"
        shell:
            """
            quast {input.polished_flye_H} {input.polished_hSpades_H} -o {params.outdir} -t {threads}
            mv {params.outdir}/report.html {output.html}
            mv {params.outdir}/report.tsv {output.tsv}
            """
