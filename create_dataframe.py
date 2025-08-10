import pathlib
import yaml
import pandas as pd
import numpy as np
from glob import glob
import os
import sys

# Load configuration from YAML file passed as the first command-line argument
with open(sys.argv[1]) as f_:
    config = yaml.load(f_, Loader=yaml.FullLoader)


def create_output_dir(output_dir):
    """
    Ensure that the output directory is located within 'results_metareads'.
    If the given path does not start with this base directory, prepend it.
    Creates the directory if it does not exist.
    
    Args:
        output_dir (str): Desired output directory (relative or full).
    
    Returns:
        str: Normalized output directory path inside 'results_metareads'.
    """
    base_dir = "results_metareads"
    if not output_dir.startswith(base_dir):
        output_dir = os.path.join(base_dir, output_dir)
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    return output_dir


def create_dataframe_illumina(input_dir):
    """
    Create a DataFrame for Illumina paired-end sequencing data.
    Expects files named with *_R1.fastq.gz and *_R2.fastq.gz.

    Args:
        input_dir (str): Path to directory containing FASTQ files.

    Returns:
        pandas.DataFrame: DataFrame with sample names and paths to R1, R2 files.
    """
    files_r1 = sorted(glob(os.path.join(input_dir, "*_R1.fastq.gz")))
    files_r2 = sorted(glob(os.path.join(input_dir, "*_R2.fastq.gz")))

    # Extract sample names from R1 and R2 filenames
    samples_r1 = {os.path.basename(f).split('_R1')[0] for f in files_r1}
    samples_r2 = {os.path.basename(f).split('_R2')[0] for f in files_r2}
    all_samples = samples_r1 | samples_r2

    # Initialize file mapping per sample
    file_mapping = {sample: {} for sample in all_samples}
    for f in files_r1:
        sample = os.path.basename(f).split('_R1')[0]
        file_mapping[sample]['_R1'] = f
    for f in files_r2:
        sample = os.path.basename(f).split('_R2')[0]
        file_mapping[sample]['_R2'] = f

    # Check for missing pairs
    missing_files = []
    for sample, files in file_mapping.items():
        for ftype in ["_R1", "_R2"]:
            if ftype not in files or not os.path.exists(files[ftype]):
                missing_files.append(f"{sample} ({ftype})")

    if missing_files:
        print("Error: The following files are missing for the illumina data type:")
        print("Existing files:")
        for missing in missing_files:
            print(f"  - {missing}")
        sys.exit(1)

    # Construct DataFrame
    data = []
    for sample, files in file_mapping.items():
        row = {
            "sample": sample,
            "fq1": files.get('_R1', np.nan),
            "fq2": files.get('_R2', np.nan),
            "ONT": np.nan
        }
        data.append(row)

    return pd.DataFrame(data)


def create_dataframe_hybrid(input_dir):
    """
    Create a DataFrame for hybrid sequencing data (Illumina + Nanopore).
    Expects *_R1.fastq.gz, *_R2.fastq.gz and *_ONT.fastq.gz files.

    Args:
        input_dir (str): Path to directory containing FASTQ files.

    Returns:
        pandas.DataFrame: DataFrame with sample names and paths to all file types.
    """
    files_r1 = sorted(glob(os.path.join(input_dir, "*_R1.fastq.gz")))
    files_r2 = sorted(glob(os.path.join(input_dir, "*_R2.fastq.gz")))
    files_ont = sorted(glob(os.path.join(input_dir, "*_ONT.fastq.gz")))

    # Extract sample names
    samples_r1 = {os.path.basename(f).split('_R1')[0] for f in files_r1}
    samples_r2 = {os.path.basename(f).split('_R2')[0] for f in files_r2}
    samples_ont = {os.path.basename(f).split('_ONT')[0] for f in files_ont}
    all_samples = samples_r1 | samples_r2 | samples_ont

    # Map file paths by sample
    file_mapping = {sample: {} for sample in all_samples}
    for f in files_r1:
        sample = os.path.basename(f).split('_R1')[0]
        file_mapping[sample]['_R1'] = f
    for f in files_r2:
        sample = os.path.basename(f).split('_R2')[0]
        file_mapping[sample]['_R2'] = f
    for f in files_ont:
        sample = os.path.basename(f).split('_ONT')[0]
        file_mapping[sample]['_ONT'] = f

    # Check for missing file types
    missing_files = []
    for sample, files in file_mapping.items():
        for ftype in ["_R1", "_R2", "_ONT"]:
            if ftype not in files or not os.path.exists(files[ftype]):
                missing_files.append(f"{sample} ({ftype})")

    if missing_files:
        print("Error: The following files are missing for the hybrid data type:")
        print("Existing files:")
        for missing in missing_files:
            print(f"  - {missing}")
        sys.exit(1)

    # Build DataFrame rows
    data = []
    for sample, files in file_mapping.items():
        row = {
            "sample": sample,
            "fq1": files.get('_R1', np.nan),
            "fq2": files.get('_R2', np.nan),
            "ONT": files.get('_ONT', np.nan)
        }
        data.append(row)

    return pd.DataFrame(data)


def create_dataframe_nanopore(input_dir):
    """
    Create a DataFrame for Nanopore-only data.
    Expects *_ONT.fastq.gz files only.

    Args:
        input_dir (str): Path to directory containing Nanopore files.

    Returns:
        pandas.DataFrame: DataFrame with sample names and ONT paths.
    """
    files_ont = sorted(glob(os.path.join(input_dir, "*_ONT.fastq.gz")))

    if not files_ont:
        print("Error: No _ONT files found. Please check the filenames in the folder.")
        sys.exit(1)

    data = []
    for f in files_ont:
        sample = os.path.basename(f).split('_ONT')[0]
        data.append({"sample": sample, "fq1": np.nan, "fq2": np.nan, "ONT": f})

    return pd.DataFrame(data)


if __name__ == '__main__':
    # Read general settings and file paths from config
    input_dir = config['general']['filename']
    data_type = config['settings']['data_type']
    output_dir = config["general"]["output_dir"]

    # Ensure output directory exists in correct base folder
    output_dir = create_output_dir(output_dir)

    # Choose appropriate DataFrame constructor based on data type
    if data_type == "illumina":
        df = create_dataframe_illumina(input_dir)
    elif data_type == "hybrid":
        df = create_dataframe_hybrid(input_dir)
    elif data_type == "nanopore":
        df = create_dataframe_nanopore(input_dir)
    else:
        print(f"Error: Unsupported data type '{data_type}'.")
        sys.exit(1)

    # Save the resulting DataFrame to a TSV file
    output_file = os.path.join(output_dir, config["general"]['units'])
    df.to_csv(output_file, sep='\t', index=False)

    # Utility function to shorten file paths for display
    def shorten_path(path, levels=2):
        """
        Return only the last 'levels' parts of a file path.
        Used to simplify output display in the terminal.
        """
        if pd.notnull(path):
            return os.path.join(*path.split(os.sep)[-levels:])
        return path

    # Create a display version of the DataFrame with shortened paths
    df_display = df.copy()
    for col in ['fq1', 'fq2', 'ONT']:
        if col in df_display.columns:
            df_display[col] = df_display[col].apply(lambda x: shorten_path(x, levels=2))

    # Print readable summary table to terminal
    print("\nmetareads; bioinformatics")
    print(df_display.dropna(how="all", axis=1))
    print(f"\nmetareads: Dataframe saved to {output_file}\n")
