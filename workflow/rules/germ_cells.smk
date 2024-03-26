rule get_adult_germ_cells:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        rds = "results/subsets/adult.sce.integrated.clustered.celltypes.germ_cell.rds",
    script:
        "../scripts/subsets/get-adult-germ-cells.R"