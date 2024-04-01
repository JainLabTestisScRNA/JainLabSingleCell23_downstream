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

rule plot_cluster_label_correlation_matrices:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
    output:
        jain_vs_jain = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.jain_vs_jain.pdf',
        jain_vs_green2018 = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.jain_vs_green2018.pdf',
        green2018_vs_green2018 = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.green2018_vs_green2018.pdf',
        no_genotype = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.no_genotype.tsv.gz',
        mut_vs_mut = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.mut_vs_mut.pdf',
        mut_vs_wt = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.mut_vs_wt.pdf',
        wt_vs_wt = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.wt_vs_wt.pdf',
        genotype = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.genotype.tsv.gz'
    params:
        odir = 'results/plots/cluster_label_correlation_matrices'
    script:
        "../scripts/plots/cluster_label_correlation_matrices.R"