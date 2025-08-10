# Rule: [dag_plot_IO]; for all data_types
rule create_dag_plot:
    priority: 1000
    output:
        report(os.path.join("results_metareads", config["general"]["output_dir"], "dag_plot_tree.pdf"), category="DAGPlot")
    run:
        dag_output = output[0]
        config_path = config["general"]["dag_config"]
        shell("snakemake --configfile {config_path} --dag | dot -Tpdf -o {dag_output}")
        print("metareads: Build DAG visualization based on config option...")
        print(f"metareads: Save DAG plot to: {dag_output}")
