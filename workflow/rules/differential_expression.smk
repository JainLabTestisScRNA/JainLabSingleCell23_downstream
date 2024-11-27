rule adult_cmo_de:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds
    output:
        pseudobulk = "results/differential_expression/adult.sce.pseudobulk.rds",
        tsv = "results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz",
        de_results = "results/differential_expression/adult.df.pseudobulk.de_results.rds"
    script:
        "../scripts/differential_expression/adult_cmo_de.R"

rule adult_cmo_germ_cell_de:
    input:
        sce = rules.lift_over_mut_adult_germ_cells.output.rds
    output:
        pseudobulk = "results/differential_expression/adult.sce.germ_cell.reprocessed.pseudobulk.rds",
        tsv = "results/differential_expression/adult.tbl.germ_cell.reprocessed.de.tsv.gz",
        de_results = "results/differential_expression/adult.df.germ_cell.reprocessed.de_results.rds"
    script:
        "../scripts/differential_expression/adult_cmo_germ_cell_de.R"