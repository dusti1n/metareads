# assembler.py
# This script is used to validate selected assemblers based on the data type defined in the config file.
import yaml
import sys

# Mapping of internal assembler keys to properly formatted names for display
FORMATTED_ASSEMBLERS = {
    "flye": "Flye",
    "canu": "Canu",
    "spades": "SPAdes",
    "megahit": "MEGAHIT"
}

# Format assembler name for output (fallback: capitalize unknown assembler)
def format_assembler_name(assembler):
    return FORMATTED_ASSEMBLERS.get(assembler, assembler.capitalize())

# Validate the assemblers selected in the configuration
def validate_assembler(config):
    # Allowed assemblers per data type
    assembler_options = {
        "illumina": {"spades", "megahit"},
        "nanopore": {"flye", "canu"},
        "hybrid": {"flye", "spades"}
    }

    # Check if required config keys are present
    if "settings" not in config or "data_type" not in config["settings"]:
        print("ERROR: 'data_type' missing in config.")
        sys.exit(1)

    selected_data_type = config["settings"]["data_type"]
    
    # Check if the selected data type is valid
    if selected_data_type not in assembler_options:
        print(f"ERROR: '{selected_data_type}' is not valid.")
        sys.exit(1)

    valid_assemblers = assembler_options[selected_data_type]
    
    # Get selected assemblers where "status" is set to true
    selected_assemblers = [
        asm for asm in config.get("assembler", {}) if config["assembler"][asm].get("status", False)
    ]

    # Validate the number of selected assemblers
    if len(selected_assemblers) < 1:
        print(f"ERROR: At least 1 assembler required for data type '{selected_data_type}'.")
        sys.exit(1)
    if len(selected_assemblers) > 2:
        print(f"ERROR: A maximum of 2 assemblers are allowed for data type '{selected_data_type}'.")
        sys.exit(1)

    # Check for assemblers not allowed for the selected data type
    invalid_assemblers = [asm for asm in selected_assemblers if asm not in valid_assemblers]
    if invalid_assemblers:
        print(f"ERROR: Invalid assemblers for {selected_data_type}: {', '.join(invalid_assemblers)}")
        sys.exit(1)
