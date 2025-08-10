**metareads: Analysis Report**

This report presents the results of a bioinformatics workflow for sequence data analysis.  
The pipeline integrates quality control, assembly, and evaluation steps to process high-throughput sequencing reads efficiently and reproducibly.

If you use this workflow in your research or publication,  
please consider citing it to acknowledge the work and support further development.  

**Repository:** https://github.com/dusti1n/metareads

----

**Excerpt from config file:**

- **Input_file:** {{ snakemake.config["general"]["filename"] }}  
- **Project_name:** {{ snakemake.config["general"]["output_dir"] }}  
- **Data_type:** {{ snakemake.config["settings"]["data_type"] }}  
- **Threads:** {{ snakemake.config["settings"]["threads"] }}  
- **Memory:** {{ snakemake.config["settings"]["memory"] }}

----

**Tools used by sequencing data type:**

- **Illumina**
  - FastQC – quality control of short reads
  - MultiQC – aggregation of FastQC results
  - Cutadapt – adapter trimming
  - SPAdes / MEGAHIT – genome or metagenome assembly
  - BWA + Pilon – assembly polishing
  - QUAST – assembly quality metrics
  - BUSCO – completeness estimation

- **Nanopore**
  - NanoPlot – visualization of long-read quality
  - Porechop – adapter trimming
  - Flye / Canu – long-read assembly
  - Racon + Medaka – polishing of long-read assemblies
  - QUAST – assembly quality metrics
  - BUSCO – completeness estimation

- **Hybrid (Illumina + Nanopore)**
  - FastQC + MultiQC + NanoPlot – dual quality assessment
  - Cutadapt + Porechop – preprocessing of both read types
  - HybridSPAdes / Flye – hybrid assembly strategies
  - Racon + Medaka + Pilon – multi-stage polishing
  - QUAST – assembly quality metrics
  - BUSCO – completeness estimation
