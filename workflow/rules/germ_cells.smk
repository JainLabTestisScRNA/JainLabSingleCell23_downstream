rule get_adult_germ_cells:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        rds = "results/germ_cells/adult.sce.germ_cell.rds",
        mut_rds = "results/germ_cells/adult.sce.germ_cell.mut.rds",
        wt_rds = "results/germ_cells/adult.sce.germ_cell.wt.rds",
        wt_spg_rds = "results/germ_cells/adult.sce.germ_cell.wt.Spermatogonia.rds",
        wt_spc_rds = "results/germ_cells/adult.sce.germ_cell.wt.Spermatocyte.rds",
        wt_rspt_rds = "results/germ_cells/adult.sce.germ_cell.wt.RoundSpermatid.rds",
        wt_espt_rds = "results/germ_cells/adult.sce.germ_cell.wt.Elongating.rds",
    script:
        "../scripts/germ_cells/get-adult-germ-cells.R"

rule subcluster_wt_adult_germ_cells:
    input:
        sce = "results/germ_cells/adult.sce.germ_cell.wt.{gc_type}.rds",
    output:
        rds = "results/germ_cells/adult.sce.germ_cell.wt.{gc_type}.subclustered.rds",
    params:
        k =  lambda wc: config.get("gc_subclusters_k").get(wc.gc_type),
    script:
        "../scripts/germ_cells/reprocess-adult-germ-cells_wt.R"

rule reintegrate_wt_adult_germ_cells:
    """
    reintegrate after subclustering broad germ cell types
    """
    input:
        sces = expand("results/germ_cells/adult.sce.germ_cell.wt.{g}.subclustered.rds",g=["Spermatogonia", "Spermatocyte", "RoundSpermatid", "Elongating"]),
    output:
        rds = "results/germ_cells/adult.sce.germ_cell.wt.subclustered.reintegrated.rds",
    script:
        "../scripts/germ_cells/reintegrate-adult-germ-cells_wt.R"

rule lift_over_mut_adult_germ_cells:
    input:
        wt = "results/germ_cells/adult.sce.germ_cell.wt.subclustered.reintegrated.rds",
        mut = "results/germ_cells/adult.sce.germ_cell.mut.rds",
    output:
        rds = "results/germ_cells/adult.sce.germ_cell.both_genotypes.subclustered.reintegrated.rds",
    script:
        "../scripts/germ_cells/lift-over-mut-adult-germ-cells.R"
