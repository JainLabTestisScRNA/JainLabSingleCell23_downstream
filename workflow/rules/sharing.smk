rule compute_sharing_score:
    input:
        sce = rules.lift_over_mut_adult_germ_cells.output.rds,
    output:
        rds = 'results/sharing/sharing_score.rds',
    script:
        "../scripts/sharing/compute_sharing_score.R"

        