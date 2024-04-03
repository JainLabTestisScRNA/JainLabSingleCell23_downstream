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

rule get_adult_germ_cells_cell_cycle:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
    output:
        rds = "results/germ_cells/adult.sce.integrated.clustered.celltypes.germ_cell.preprocessed.cell_cycle.rds",
        cc = "results/germ_cells/adult.integrated.clustered.celltypes.germ_cell.cell_cycle_classification.rds",
    threads: 6,
    script:
        "../scripts/germ_cells/get-adult-germ-cells-cell-cycle.R"