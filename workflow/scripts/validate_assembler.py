import yaml
import sys

FORMATTED_ASSEMBLERS = {
    "flye": "Flye",
    "canu": "Canu",
    "spades": "SPAdes",
    "megahit": "MEGAHIT"
}

def format_assembler_name(assembler):
    return FORMATTED_ASSEMBLERS.get(assembler, assembler.capitalize())

def validate_assembler(config):
    assembler_options = {
        "Illumina": {"spades", "megahit"},
        "Nanopore": {"flye", "canu"},
        "Hybrid": {"flye", "spades"}
    }

    if "settings" not in config or "data_type" not in config["settings"]:
        print("ERROR: 'data_type' missing in config.")
        sys.exit(1)

    selected_data_type = config["settings"]["data_type"]
    if selected_data_type not in assembler_options:
        print(f"ERROR: '{selected_data_type}' is not valid.")
        sys.exit(1)

    valid_assemblers = assembler_options[selected_data_type]
    selected_assemblers = [
        asm for asm in config.get("assembler", {}) if config["assembler"][asm].get("status", False)
    ]

    if len(selected_assemblers) < 1:
        print(f"ERROR: At least 1 assembler required for data type '{selected_data_type}'.")
        sys.exit(1)
    if len(selected_assemblers) > 2:
        print(f"ERROR: A maximum of 2 assemblers are allowed for data type '{selected_data_type}'.")
        sys.exit(1)

    invalid_assemblers = [asm for asm in selected_assemblers if asm not in valid_assemblers]
    if invalid_assemblers:
        print(f"ERROR: Invalid assemblers for {selected_data_type}: {', '.join(invalid_assemblers)}")
        sys.exit(1)

