#### metareads; Bioinformatics
# Changelog
This document records all notable changes to metareads project.

## [2025-06-03]  
- Create config_presets folder  
- Add configuration files to /config_presets  
- Add DAG plot option to config_file 
- Update all included rules for thread management  
- Update and fix bugs in validate_assemlber.py

## [2025-05-13]  
- Add Illumina assembler: MEGAHIT  
- Add env file for MEGAHIT assembler  
- Update all rules and dependencies  
- Update Snakefile and dataframe.py  
- Update CHANGELOG.md

## [2025-02-19]  
- Downgrade Canu version from 2.3 to 2.2 due to Java bug 

## [2025-02-18]  
- Integrate SPAdes and HybridSPAdes assemblers  
- Update config.yaml with new settings  
- Modify and add new comments  

## [2025-02-17]  
- Update config.yaml with new settings  
- Add validate_config.py to validate config.yaml

## [2025-02-14]  
- Update dataframe.py  
- Define assembler settings in config.yaml  
- Integrate Flye and Canu assemblers for ONT data  
- Rename files for consistency  
