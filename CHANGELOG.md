#### metareads; bioinformatics

# Changelog

This document records all notable changes to the **metareads** project.  
**Note**: Use the **dev** branch for the latest updates and features.

---

### [2025-10-08]
- Add post-polish workflow for all data types  
  - Tools: Racon, Pilon, Medaka
- Add post-analysis workflow for metagenomic results  
  - Tools: QUAST, BUSCO, custom R script for BUSCO plot
- Create full merge process for all data types  
  - Integrate Snakemake report interface  
  - Merge result files (e.g., .txt, .csv)  
  - Configure Snakemake report for evaluation output  
  - Add analysis summary section to the report
- Improve codebase  
  - Clean code and add meaningful comments in various files  
  - Revise root directory structure
- Rename build_dataframe.py to create_dataframe.py
- Rename assembler rule files for clarity
- Optimize assemblers for metagenomic datasets
- Update dataframe.py with function to better summarize input table
- Revise environment files
- Add preset configuration files for example datasets
- Modularize entire workflow for all data type-specific paths
- Include main metareads report site to provide metadata overview
- Add R scripts: buscoplot.R, buscosum.R
  - buscoplot.R: Generate plot after BUSCO to evaluate final assemblies  
  - buscosum.R: Combine final plots into one for comparison
- Perform test run of metareads using example data
- Update main metareads.yaml environment configuration
- Create README.md and add usage guide for metareads
- Mark selected output files as temp() to reduce disk usage

### [2025-06-03]  
- Create config_presets folder  
- Add configuration files to /config_presets
- Add DAG plot option to config_file 
- Update all included rules for thread management  
- Fix bugs and update validate_assembler.py

### [2025-05-13]  
- Add Illumina assembler: MEGAHIT  
- Add environment file for MEGAHIT assembler  
- Update all rules and dependencies  
- Update Snakefile and dataframe.py
- Update CHANGELOG.md

### [2025-02-19]  
- Downgrade Canu version from 2.3 to 2.2 due to Java bug

### [2025-02-18]  
- Integrate SPAdes and HybridSPAdes assemblers  
- Update config.yaml with new settings  
- Modify and add new comments

### [2025-02-17]  
- Update config.yaml with new settings  
- Add validate_config.py to validate config.yaml

### [2025-02-14]  
- Update dataframe.py
- Define assembler settings in config.yaml
- Integrate Flye and Canu assemblers for ONT data  
- Rename files for consistency
