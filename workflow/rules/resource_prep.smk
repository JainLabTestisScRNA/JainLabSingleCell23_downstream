rule get_cell_cycle_genes:
    input:
        xlsx=config.get("cell_cycle_genes")
    output:
        tsv="results/resources/cell_cycle_genes.tsv"
    script:
        "../scripts/resource_prep/cell_cycle_genes.R"