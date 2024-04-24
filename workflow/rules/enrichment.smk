rule gsea_prerank:
    input:
        sce = rules.lift_over_mut_adult_germ_cells.output.rds,
        msigdb = config.get("msigdb"),
    output:
        rds = "results/enrichment/gsea_prerank.rds",
    script:
        "../scripts/enrichment/gsea_prerank.R"