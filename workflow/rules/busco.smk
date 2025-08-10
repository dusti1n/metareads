# Data type: Illumina
# Pathway (Spades/MEGAHIT): Spades/MEGAHIT-BWA-Pilon-(BUSCO)
# Analyze final assemblies with BUSCO

# SECTION: COMPLETE BUSCO ANALYSIS FOR ILLUMINA DATASET
if config["settings"]["data_type"] == "illumina":
# SECTION_1_ILLUMINA: BUSCO ANALYSIS
    # Pathway: Spades-BWA-Pilon-(BUSCO)
    # Rule: [busco_analysis_spades_I]; for illumina data type
    rule busco_analysis_spades_I:
        input:
            # Pathway (Spades): Spades-BWA-Pilon-(BUSCO)
            # Input_file: {sample}_pilon_spades_I_final_assemb.fasta
            polished_spades_I="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_pilon_spades_I_final_assemb.fasta"
        output:
            txt=report("results_metareads/{project_name}/analysis/busco/spades/{sample}_busco_sum_spades_I.txt", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        threads: config["settings"]["threads"]
        params:
            tmp_dir="{sample}_spades_I",
            mode="genome",
            lineage=config["settings"]["busco_lineage"],
            output_path="results_metareads/{project_name}/analysis/busco/spades",
            download_path="results_metareads/{project_name}/analysis/busco/lineage_downloads"
        shell:
            """
            busco -i {input.polished_spades_I} \
                -o {params.tmp_dir} \
                -l {params.lineage} \
                -m {params.mode} \
                -c  {threads} \
                --out_path {params.output_path} \
                --download_path {params.download_path}

            cp {params.output_path}/{params.tmp_dir}/short_summary.specific.{params.lineage}.{params.tmp_dir}.txt \
            {output.txt}
            """

    # Pathway: MEGAHIT-BWA-Pilon-(BUSCO)
    # Rule: [busco_analysis_megahit_I]; for illumina data type
    rule busco_analysis_megahit_I:
        input:
            # Pathway (MEGAHIT): MEGAHIT-BWA-Pilon-(BUSCO)
            # Input_file: {sample}_pilon_megahit_I_final_assemb.fasta
            polished_megahit_I="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_pilon_megahit_I_final_assemb.fasta"
        output:
            txt=report("results_metareads/{project_name}/analysis/busco/megahit/{sample}_busco_sum_megahit_I.txt", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        threads: config["settings"]["threads"]
        params:
            tmp_dir="{sample}_megahit_I",
            mode="genome",
            lineage=config["settings"]["busco_lineage"],
            output_path="results_metareads/{project_name}/analysis/busco/megahit",
            download_path="results_metareads/{project_name}/analysis/busco/lineage_downloads"
        shell:
            """
            busco -i {input.polished_megahit_I} \
                -o {params.tmp_dir} \
                -l {params.lineage} \
                -m {params.mode} \
                -c  {threads} \
                --out_path {params.output_path} \
                --download_path {params.download_path}

            cp {params.output_path}/{params.tmp_dir}/short_summary.specific.{params.lineage}.{params.tmp_dir}.txt \
            {output.txt}
            """


# SECTION_2_ILLUMINA: BUSCO WRITE CSV
    # Pathway: Spades-BWA-Pilon-(BUSCO)
    # Rule: [busco_csv_spades_I]; for illumina data type
    rule busco_csv_spades_I:
            input:
                busco_sum="results_metareads/{project_name}/analysis/busco/spades/{sample}_busco_sum_spades_I.txt"
            output:
                csv_file="results_metareads/{project_name}/analysis/busco/spades/{sample}_busco_sum_spades_I.csv"
            conda:
                "../envs/busco.yaml"
            threads: config["settings"]["threads"]
            params:
                temp_dir="results_metareads/{project_name}/analysis/busco/spades/{sample}_temp/",
                prefix="{sample}_spades_I"
            shell:
                """
                echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.csv_file}
                mkdir -p {params.temp_dir}
                cat {input.busco_sum} | grep "(S)" | awk -v strain="{params.prefix}" '{{print strain","$1}}' > {params.temp_dir}/complete_single.txt
                cat {input.busco_sum} | grep "(D)" | awk '{{print $1}}' > {params.temp_dir}/complete_duplicated.txt
                cat {input.busco_sum} | grep "(F)" | awk '{{print $1}}' > {params.temp_dir}/fragmented.txt
                cat {input.busco_sum} | grep "(M)" | awk '{{print $1}}' > {params.temp_dir}/missing.txt
                paste -d "," {params.temp_dir}/complete_single.txt {params.temp_dir}/complete_duplicated.txt {params.temp_dir}/fragmented.txt {params.temp_dir}/missing.txt >> {output.csv_file}
                rm -r {params.temp_dir}
                """

    # Pathway: MEGAHIT-BWA-Pilon-(BUSCO)
    # Rule: [busco_csv_megahit_I]; for illumina data type
    rule busco_csv_megahit_I:
            input:
                busco_sum="results_metareads/{project_name}/analysis/busco/megahit/{sample}_busco_sum_megahit_I.txt"
            output:
                csv_file="results_metareads/{project_name}/analysis/busco/megahit/{sample}_busco_sum_megahit_I.csv"
            conda:
                "../envs/busco.yaml"
            threads: config["settings"]["threads"]
            params:
                temp_dir="results_metareads/{project_name}/analysis/busco/megahit/{sample}_temp/",
                prefix="{sample}_megahit_I"
            shell:
                """
                echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.csv_file}
                mkdir -p {params.temp_dir}
                cat {input.busco_sum} | grep "(S)" | awk -v strain="{params.prefix}" '{{print strain","$1}}' > {params.temp_dir}/complete_single.txt
                cat {input.busco_sum} | grep "(D)" | awk '{{print $1}}' > {params.temp_dir}/complete_duplicated.txt
                cat {input.busco_sum} | grep "(F)" | awk '{{print $1}}' > {params.temp_dir}/fragmented.txt
                cat {input.busco_sum} | grep "(M)" | awk '{{print $1}}' > {params.temp_dir}/missing.txt
                paste -d "," {params.temp_dir}/complete_single.txt {params.temp_dir}/complete_duplicated.txt {params.temp_dir}/fragmented.txt {params.temp_dir}/missing.txt >> {output.csv_file}
                rm -r {params.temp_dir}
                """


# SECTION_3_ILLUMINA: BUSCO CREATE PLOT
    # Pathway: Spades-BWA-Pilon-(BUSCO)
    # Rule: [busco_plot_spades_I]; for illumina data type
    rule busco_plot_spades_I:
        input:
            csv_file="results_metareads/{project_name}/analysis/busco/spades/{sample}_busco_sum_spades_I.csv"
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/spades/{sample}_busco_plot_spades_I.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscoplot.R {input.csv_file} {output.plot} {wildcards.sample}
            """

    # Pathway: MEGAHIT-BWA-Pilon-(BUSCO)
    # Rule: [busco_plot_megahit_I]; for illumina data type
    rule busco_plot_megahit_I:
        input:
            csv_file="results_metareads/{project_name}/analysis/busco/megahit/{sample}_busco_sum_megahit_I.csv"
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/megahit/{sample}_busco_plot_megahit_I.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscoplot.R {input.csv_file} {output.plot} {wildcards.sample}
            """


# SECTION_4_ILLUMINA: BUSCO PLOT SUM
    # Pathway: Spades/MEGAHIT-BWA-Pilon-(BUSCO)
    # Rule: [busco_plot_sum_I]; for illumina data type
    rule busco_plot_sum_I:
        input:
            # Pathway (Spades): Spades-BWA-Pilon-(BUSCO)
            # Input_file: {sample}_busco_plot_spades_I.pdf
            busco_plot_spades_I="results_metareads/{project_name}/analysis/busco/spades/{sample}_busco_sum_spades_I.csv"
                if config["settings"]["data_type"] == "illumina" and config["assembler"]["spades"]["status"] else [],

            # Pathway (MEGAHIT): MEGAHIT-BWA-Pilon-(BUSCO)
            # Input_file: {sample}_busco_plot_megahit_I.pdf
            busco_plot_megahit_I="results_metareads/{project_name}/analysis/busco/megahit/{sample}_busco_sum_megahit_I.csv"
                if config["settings"]["data_type"] == "illumina" and config["assembler"]["megahit"]["status"] else [],
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/sum/{sample}_busco_plot_sum_I.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscosum.R {input} {output.plot}
            """




# Data type: Nanopore
# Pathway (Flye/Canu): Flye/Canu-Minimap2-Racon-MEDAKA-(BUSCO)
# Analyze final assemblies with BUSCO

# SECTION: COMPLETE BUSCO ANALYSIS FOR NANOPORE DATASET
if config["settings"]["data_type"] == "nanopore":
# SECTION_1_NANOPORE: BUSCO ANALYSIS
    # Pathway: Flye-Minimap2-Racon-MEDAKA-(BUSCO)
    # Rule: [busco_analysis_flye_N]; for nanopore data type
    rule busco_analysis_flye_N:
        input:
            # Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-(BUSCO)
            # Input_file: {sample}_flye_N_final_assemb.fasta
            polished_flye_medaka_N="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_flye_N_final_assemb.fasta"
        output:
            txt=report("results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_N.txt", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        threads: config["settings"]["threads"]
        params:
            tmp_dir="{sample}_flye_N",
            mode="genome",
            lineage=config["settings"]["busco_lineage"],
            output_path="results_metareads/{project_name}/analysis/busco/flye",
            download_path="results_metareads/{project_name}/analysis/busco/lineage_downloads"
        shell:
            """
            busco -i {input.polished_flye_medaka_N} \
                -o {params.tmp_dir} \
                -l {params.lineage} \
                -m {params.mode} \
                -c  {threads} \
                --out_path {params.output_path} \
                --download_path {params.download_path}

            cp {params.output_path}/{params.tmp_dir}/short_summary.specific.{params.lineage}.{params.tmp_dir}.txt \
            {output.txt}
            """

    # Pathway: Canu-Minimap2-Racon-MEDAKA-(BUSCO)
    # Rule: [busco_analysis_canu_N]; for nanopore data type
    rule busco_analysis_canu_N:
        input:
            # Pathway (Canu): Canu-Minimap2-Racon-MEDAKA-(BUSCO)
            # Input_file: {sample}_canu_N_final_assemb.fasta
            polished_canu_medaka_N="results_metareads/{project_name}/filtered/medaka/{sample}/{sample}_canu_N_final_assemb.fasta"
        output:
            txt=report("results_metareads/{project_name}/analysis/busco/canu/{sample}_busco_sum_canu_N.txt", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        threads: config["settings"]["threads"]
        params:
            tmp_dir="{sample}_canu_N",
            mode="genome",
            lineage=config["settings"]["busco_lineage"],
            output_path="results_metareads/{project_name}/analysis/busco/canu",
            download_path="results_metareads/{project_name}/analysis/busco/lineage_downloads"
        shell:
            """
            busco -i {input.polished_canu_medaka_N} \
                -o {params.tmp_dir} \
                -l {params.lineage} \
                -m {params.mode} \
                -c  {threads} \
                --out_path {params.output_path} \
                --download_path {params.download_path}

            cp {params.output_path}/{params.tmp_dir}/short_summary.specific.{params.lineage}.{params.tmp_dir}.txt \
            {output.txt}
            """


# SECTION_2_NANOPORE: BUSCO WRITE CSV
    # Pathway: Flye-Minimap2-Racon-MEDAKA-(BUSCO)
    # Rule: [busco_csv_flye_N]; for nanopore data type
    rule busco_csv_flye_N:
            input:
                busco_sum="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_N.txt"
            output:
                csv_file="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_N.csv"
            conda:
                "../envs/busco.yaml"
            threads: config["settings"]["threads"]
            params:
                temp_dir="results_metareads/{project_name}/analysis/busco/flye/{sample}_temp/",
                prefix="{sample}_flye_N"
            shell:
                """
                echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.csv_file}
                mkdir -p {params.temp_dir}
                cat {input.busco_sum} | grep "(S)" | awk -v strain="{params.prefix}" '{{print strain","$1}}' > {params.temp_dir}/complete_single.txt
                cat {input.busco_sum} | grep "(D)" | awk '{{print $1}}' > {params.temp_dir}/complete_duplicated.txt
                cat {input.busco_sum} | grep "(F)" | awk '{{print $1}}' > {params.temp_dir}/fragmented.txt
                cat {input.busco_sum} | grep "(M)" | awk '{{print $1}}' > {params.temp_dir}/missing.txt
                paste -d "," {params.temp_dir}/complete_single.txt {params.temp_dir}/complete_duplicated.txt {params.temp_dir}/fragmented.txt {params.temp_dir}/missing.txt >> {output.csv_file}
                rm -r {params.temp_dir}
                """

    # Pathway: Canu-Minimap2-Racon-MEDAKA-(BUSCO)
    # Rule: [busco_csv_canu_N]; for nanopore data type
    rule busco_csv_canu_N:
            input:
                busco_sum="results_metareads/{project_name}/analysis/busco/canu/{sample}_busco_sum_canu_N.txt"
            output:
                csv_file="results_metareads/{project_name}/analysis/busco/canu/{sample}_busco_sum_canu_N.csv"
            conda:
                "../envs/busco.yaml"
            threads: config["settings"]["threads"]
            params:
                temp_dir="results_metareads/{project_name}/analysis/busco/canu/{sample}_temp/",
                prefix="{sample}_canu_N"
            shell:
                """
                echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.csv_file}
                mkdir -p {params.temp_dir}
                cat {input.busco_sum} | grep "(S)" | awk -v strain="{params.prefix}" '{{print strain","$1}}' > {params.temp_dir}/complete_single.txt
                cat {input.busco_sum} | grep "(D)" | awk '{{print $1}}' > {params.temp_dir}/complete_duplicated.txt
                cat {input.busco_sum} | grep "(F)" | awk '{{print $1}}' > {params.temp_dir}/fragmented.txt
                cat {input.busco_sum} | grep "(M)" | awk '{{print $1}}' > {params.temp_dir}/missing.txt
                paste -d "," {params.temp_dir}/complete_single.txt {params.temp_dir}/complete_duplicated.txt {params.temp_dir}/fragmented.txt {params.temp_dir}/missing.txt >> {output.csv_file}
                rm -r {params.temp_dir}
                """


# SECTION_3_NANOPORE: BUSCO CREATE PLOT
    # Pathway: Flye-Minimap2-Racon-MEDAKA-(BUSCO)
    # Rule: [busco_plot_flye_N]; for nanopore data type
    rule busco_plot_flye_N:
        input:
            csv_file="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_N.csv"
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_plot_flye_N.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscoplot.R {input.csv_file} {output.plot} {wildcards.sample}
            """

    # Pathway: Canu-Minimap2-Racon-MEDAKA-(BUSCO)
    # Rule: [busco_plot_canu_N]; for nanopore data type
    rule busco_plot_canu_N:
        input:
            csv_file="results_metareads/{project_name}/analysis/busco/canu/{sample}_busco_sum_canu_N.csv"
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/canu/{sample}_busco_plot_canu_N.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscoplot.R {input.csv_file} {output.plot} {wildcards.sample}
            """


# SECTION_4_NANOPORE: BUSCO PLOT SUM
    # Pathway: Flye/Canu-Minimap2-Racon-MEDAKA-(BUSCO)
    # Rule: [busco_plot_sum_N]; for nanopore data type
    rule busco_plot_sum_N:
        input:
            # Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-(BUSCO)
            # Input_file: {sample}_busco_sum_flye_N.csv
            busco_plot_flye_N="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_N.csv"
                if config["settings"]["data_type"] == "nanopore" and config["assembler"]["flye"]["status"] else [],

            # Pathway (Canu): Canu-Minimap2-Racon-MEDAKA-(BUSCO)
            # Input_file: {sample}_busco_sum_canu_N.csv
            busco_plot_canu_N="results_metareads/{project_name}/analysis/busco/canu/{sample}_busco_sum_canu_N.csv"
                if config["settings"]["data_type"] == "nanopore" and config["assembler"]["canu"]["status"] else [],
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/sum/{sample}_busco_plot_sum_N.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscosum.R {input} {output.plot}
            """




# Data type: Hybrid
# 1_Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-BWA-Pilon-(BUSCO)
# 2_Pathway (hSpades): hSpades-BWA-Pilon-(BUSCO)
# Analyze final assemblies with BUSCO

# SECTION: COMPLETE BUSCO ANALYSIS FOR HYBRID DATASET
if config["settings"]["data_type"] == "hybrid":
# SECTION_1_HYBRID: BUSCO ANALYSIS
    # 1_Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-BWA-Pilon-(BUSCO)
    # Rule: [busco_analysis_flye_H]; for hybrid data type
    rule busco_analysis_flye_H:
        input:
            # Input_file: {sample}_flye_H_final_assemb.fasta
            polished_flye_H="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_flye_H_final_assemb.fasta"
        output:
            txt=report("results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_H.txt", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        threads: config["settings"]["threads"]
        params:
            tmp_dir="{sample}_flye_H",
            mode="genome",
            lineage=config["settings"]["busco_lineage"],
            output_path="results_metareads/{project_name}/analysis/busco/flye",
            download_path="results_metareads/{project_name}/analysis/busco/lineage_downloads"
        shell:
            """
            busco -i {input.polished_flye_H} \
                -o {params.tmp_dir} \
                -l {params.lineage} \
                -m {params.mode} \
                -c  {threads} \
                --out_path {params.output_path} \
                --download_path {params.download_path}

            cp {params.output_path}/{params.tmp_dir}/short_summary.specific.{params.lineage}.{params.tmp_dir}.txt \
            {output.txt}
            """

    # 2_Pathway (hSpades): hSpades-BWA-Pilon-(BUSCO)
    # Rule: [busco_analysis_hSpades_H]; for hybrid data type
    rule busco_analysis_hSpades_H:
        input:
            # Input_file: {sample}_hSpades_H_final_assemb.fasta
            polished_hSpades_H="results_metareads/{project_name}/filtered/pilon/{sample}/{sample}_hSpades_H_final_assemb.fasta"
        output:
            txt=report("results_metareads/{project_name}/analysis/busco/hSpades/{sample}_busco_sum_hSpades_H.txt", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        threads: config["settings"]["threads"]
        params:
            tmp_dir="{sample}_hSpades_H",
            mode="genome",
            lineage=config["settings"]["busco_lineage"],
            output_path="results_metareads/{project_name}/analysis/busco/hSpades",
            download_path="results_metareads/{project_name}/analysis/busco/lineage_downloads"
        shell:
            """
            busco -i {input.polished_hSpades_H} \
                -o {params.tmp_dir} \
                -l {params.lineage} \
                -m {params.mode} \
                -c  {threads} \
                --out_path {params.output_path} \
                --download_path {params.download_path}

            cp {params.output_path}/{params.tmp_dir}/short_summary.specific.{params.lineage}.{params.tmp_dir}.txt \
            {output.txt}
            """


# SECTION_2_HYBRID: BUSCO WRITE CSV
    # 1_Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-BWA-Pilon-(BUSCO)
    # Rule: [busco_csv_flye_H]; for hybrid data type
    rule busco_csv_flye_H:
            input:
                busco_sum="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_H.txt"
            output:
                csv_file="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_H.csv"
            conda:
                "../envs/busco.yaml"
            threads: config["settings"]["threads"]
            params:
                temp_dir="results_metareads/{project_name}/analysis/busco/flye/{sample}_temp/",
                prefix="{sample}_flye_H"
            shell:
                """
                echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.csv_file}
                mkdir -p {params.temp_dir}
                cat {input.busco_sum} | grep "(S)" | awk -v strain="{params.prefix}" '{{print strain","$1}}' > {params.temp_dir}/complete_single.txt
                cat {input.busco_sum} | grep "(D)" | awk '{{print $1}}' > {params.temp_dir}/complete_duplicated.txt
                cat {input.busco_sum} | grep "(F)" | awk '{{print $1}}' > {params.temp_dir}/fragmented.txt
                cat {input.busco_sum} | grep "(M)" | awk '{{print $1}}' > {params.temp_dir}/missing.txt
                paste -d "," {params.temp_dir}/complete_single.txt {params.temp_dir}/complete_duplicated.txt {params.temp_dir}/fragmented.txt {params.temp_dir}/missing.txt >> {output.csv_file}
                rm -r {params.temp_dir}
                """

    # 2_Pathway (hSpades): hSpades-BWA-Pilon-(BUSCO)
    # Rule: [busco_csv_hSpades_H]; for hybrid data type
    rule busco_csv_hSpades_H:
            input:
                busco_sum="results_metareads/{project_name}/analysis/busco/hSpades/{sample}_busco_sum_hSpades_H.txt"
            output:
                csv_file="results_metareads/{project_name}/analysis/busco/hSpades/{sample}_busco_sum_hSpades_H.csv"
            conda:
                "../envs/busco.yaml"
            threads: config["settings"]["threads"]
            params:
                temp_dir="results_metareads/{project_name}/analysis/busco/hSpades/{sample}_temp/",
                prefix="{sample}_hSpades_H"
            shell:
                """
                echo "Strain,Complete_single_copy,Complete_duplicated,Fragmented,Missing" > {output.csv_file}
                mkdir -p {params.temp_dir}
                cat {input.busco_sum} | grep "(S)" | awk -v strain="{params.prefix}" '{{print strain","$1}}' > {params.temp_dir}/complete_single.txt
                cat {input.busco_sum} | grep "(D)" | awk '{{print $1}}' > {params.temp_dir}/complete_duplicated.txt
                cat {input.busco_sum} | grep "(F)" | awk '{{print $1}}' > {params.temp_dir}/fragmented.txt
                cat {input.busco_sum} | grep "(M)" | awk '{{print $1}}' > {params.temp_dir}/missing.txt
                paste -d "," {params.temp_dir}/complete_single.txt {params.temp_dir}/complete_duplicated.txt {params.temp_dir}/fragmented.txt {params.temp_dir}/missing.txt >> {output.csv_file}
                rm -r {params.temp_dir}
                """


# SECTION_3_HYBRID: BUSCO CREATE PLOT
    # 1_Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-BWA-Pilon-(BUSCO)
    # Rule: [busco_plot_flye_H]; for hybrid data type
    rule busco_plot_flye_H:
        input:
            csv_file="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_H.csv"
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_plot_flye_H.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscoplot.R {input.csv_file} {output.plot} {wildcards.sample}
            """

    # 2_Pathway (hSpades): hSpades-BWA-Pilon-(BUSCO)
    # Rule: [busco_plot_hSpades_H]; for hybrid data type
    rule busco_plot_hSpades_H:
        input:
            csv_file="results_metareads/{project_name}/analysis/busco/hSpades/{sample}_busco_sum_hSpades_H.csv"
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/hSpades/{sample}_busco_plot_hSpades_H.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscoplot.R {input.csv_file} {output.plot} {wildcards.sample}
            """


# SECTION_4_HYBRID: BUSCO PLOT SUM
    # 1_Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-BWA-Pilon-(BUSCO)
    # 2_Pathway (hSpades): hSpades-BWA-Pilon-(BUSCO)
    # Rule: [busco_plot_sum_H]; for hybrid data type
    rule busco_plot_sum_H:
        input:
            # 1_Pathway (Flye): Flye-Minimap2-Racon-MEDAKA-BWA-Pilon-(BUSCO)
            # Input_file: {sample}_busco_sum_flye_H.csv
            busco_plot_flye_H="results_metareads/{project_name}/analysis/busco/flye/{sample}_busco_sum_flye_H.csv"
                if config["settings"]["data_type"] == "hybrid" and config["assembler"]["flye"]["status"] else [],

            # 2_Pathway (hSpades): hSpades-BWA-Pilon-(BUSCO)
            # Input_file: {sample}_busco_sum_hSpades_H.csv
            busco_plot_hSpades_H="results_metareads/{project_name}/analysis/busco/hSpades/{sample}_busco_sum_hSpades_H.csv"
                if config["settings"]["data_type"] == "hybrid" and config["assembler"]["spades"]["status"] else [],
        output:
            plot=report("results_metareads/{project_name}/analysis/busco/sum/{sample}_busco_plot_sum_H.pdf", category="BUSCO")
        conda:
            "../envs/busco.yaml"
        shell:
            """
            Rscript workflow/scripts/buscosum.R {input} {output.plot}
            """
