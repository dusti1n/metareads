rule build_dag_plot:
    output:
        os.path.join("results", config["general"]["output_dir"], "dag_full_io.pdf")
    run:
        dag_output = output[0]
        config_path = config["general"]["dag_config"]
        shell("snakemake --configfile {config_path} --dag | dot -Tpdf -o {dag_output}")
        print("[ms] Build DAG visualization based on config option...")
        print(f"[ms] Save DAG plot to: {dag_output}")
