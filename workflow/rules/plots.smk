rule plot_germ_and_somatic_pca:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        g_pc1_pc2 = 'results/plots/germ_and_somatic_pca/germ_and_somatic_pca.pc1_pc2.pdf',
        g_pc1_pc3 = 'results/plots/germ_and_somatic_pca/germ_and_somatic_pca.pc1_pc3.pdf',
        tsv = 'results/plots/germ_and_somatic_pca/germ_and_somatic_pca.tsv.gz'
    params:
        odir = 'results/plots/germ_and_somatic_pca'
    script:
        "../scripts/plots/germ_and_somatic_pca.R"