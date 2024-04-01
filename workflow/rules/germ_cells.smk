rule get_adult_germ_cells:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        rds = "results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.rds",
    script:
        "../scripts/germ_cells/get-adult-germ-cells.R"

rule reprocess_adult_germ_cells:
    input:
        sce = rules.get_adult_germ_cells.output.rds,
    output:
        rds = "results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.reprocessed.rds",
    script:
        "../scripts/germ_cells/reprocess-adult-germ-cells.R"