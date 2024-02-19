rule adult_cmo_de:
    input:
        sce = rules.integrate_adult_cmo.output.rds
    output:
        pseudobulk = "results/differential_expression/adult.sce.pseudobulk.rds",
        tsv = "results/differential_expression/adult.tbl.pseudobulk.de.tsv.gz",
        de_results = "results/differential_expression/adult.df.pseudobulk.de_results.rds"
    script:
        "../scripts/differential_expression/adult_cmo_de.R"