rule plot_barcode_overlap:
    input:
        fls = expand("results/preprocessing/{x}.sce.passing_cells.rds",x=config.get("cellranger_matrices")),
    output:
        tsv = 'results/plots/barcode_overlap/barcode_overlap.tsv.gz',
        pdf = 'results/plots/barcode_overlap/barcode_overlap.pdf'
    script:
        "../scripts/plots/barcode_overlap.R"

rule plot_check_normalization:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        tsv = 'results/plots/check_normalization/check_normalization.tsv.gz',
        pdf = 'results/plots/check_normalization/check_normalization.pdf'
    script:
        "../scripts/plots/check_normalization.R"

rule plot_genotype_facetted_umap:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        pdf = 'results/plots/genotype_facetted_umap/genotype_facetted_umap.pdf',
        tsv = 'results/plots/genotype_facetted_umap/genotype_facetted_umap.tsv.gz'
    script:
        "../scripts/plots/genotype_facetted_umap.R"

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

rule plot_reprocessed_germ_cell_pca:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
    output:
        pc1_pc2 = 'results/plots/reprocessed_germ_cell_pca/reprocessed_germ_cell_pca.pc1_pc2.pdf',
        pc1_pc3 = 'results/plots/reprocessed_germ_cell_pca/reprocessed_germ_cell_pca.pc1_pc3.pdf',
        tsv = 'results/plots/reprocessed_germ_cell_pca/reprocessed_germ_cell_pca.tsv.gz'
    params:
        odir = 'results/plots/reprocessed_germ_cell_pca'
    script:
        "../scripts/plots/reprocessed_germ_cell_pca.R"

rule plot_cluster_label_correlation_matrices:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
    output:
        jain_vs_jain = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.jain_vs_jain.pdf',
        jain_vs_green2018 = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.jain_vs_green2018.pdf',
        mut_vs_green2018 = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.mut_vs_green2018.pdf',
        wt_vs_green2018 = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.wt_vs_green2018.pdf',
        green2018_vs_green2018 = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.green2018_vs_green2018.pdf',
        mat = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.tsv.gz',
        mut_vs_mut = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.mut_vs_mut.pdf',
        mut_vs_wt = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.mut_vs_wt.pdf',
        wt_vs_wt = 'results/plots/cluster_label_correlation_matrices/cluster_label_correlation_matrix.wt_vs_wt.pdf',
    params:
        odir = 'results/plots/cluster_label_correlation_matrices'
    script:
        "../scripts/plots/cluster_label_correlation_matrices.R"


rule plot_cluster_label_correlation_matrices_wt_mut:
    input:
        mut = rules.reprocess_adult_germ_cells_mut.output.rds,
        wt = rules.reprocess_adult_germ_cells_wt.output.rds,
    output:
        wt_vs_green2018 = 'results/plots/cluster_label_correlation_matrices_wt_mut/cluster_label_correlation_matrix.wt_vs_green2018.pdf',
        mut_vs_green2018 = 'results/plots/cluster_label_correlation_matrices_wt_mut/cluster_label_correlation_matrix.mut_vs_green2018.pdf',
        mat = 'results/plots/cluster_label_correlation_matrices_wt_mut/cluster_label_correlation_matrix.tsv.gz',
        mut_vs_wt = 'results/plots/cluster_label_correlation_matrices_wt_mut/cluster_label_correlation_matrix.mut_vs_wt.pdf',
    params:
        odir = 'results/plots/cluster_label_correlation_matrices_wt_mut'
    script:
        "../scripts/plots/cluster_label_correlation_matrices_wt_mut.R"

rule plot_overall_marker_genes:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        tsv = 'results/plots/overall_marker_genes/overall_marker_genes.tsv.gz',
        pdf = 'results/plots/overall_marker_genes/overall_marker_genes.pdf'
    script:
        "../scripts/plots/overall_marker_genes.R"

rule plot_tex14_umap:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        pdf = 'results/plots/tex14_umap/tex14_umap.pdf',
        tsv = 'results/plots/tex14_umap/tex14_umap.tsv.gz'
    script:
        "../scripts/plots/tex14_umap.R"

rule plot_tex14_violin_germ_clusters:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
    output:
        pdf = 'results/plots/tex14_violin_germ_clusters/tex14_violin_germ_clusters.pdf',
        tsv = 'results/plots/tex14_violin_germ_clusters/tex14_violin_germ_clusters.tsv.gz'
    script:
        "../scripts/plots/tex14_violin_germ_clusters.R"

rule plot_tex14_violin_all_cells:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        pdf = 'results/plots/tex14_violin_all_cells/tex14_violin_all_cells.pdf',
        tsv = 'results/plots/tex14_violin_all_cells/tex14_violin_all_cells.tsv.gz'
    script:
        "../scripts/plots/tex14_violin_all_cells.R"

rule plot_preL_violin:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
    output:
        pdf = 'results/plots/preL_violin/preL_violin.pdf',
        tsv = 'results/plots/preL_violin/preL_violin.tsv.gz'
    script:
        "../scripts/plots/preL_violin.R"

rule plot_preL_marker_pca:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
    output:
        pdf = 'results/plots/preL_marker_pca/preL_marker_pca.pdf',
        tsv = 'results/plots/preL_marker_pca/preL_marker_pca.tsv.gz'
    script:
        "../scripts/plots/preL_marker_pca.R"

rule plot_germ_cell_cell_cycle_bar:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        pdf = 'results/plots/germ_cell_cycle_bar/germ_cell_cycle_bar.pdf',
        tsv = 'results/plots/germ_cell_cycle_bar/germ_cell_cycle_bar.tsv.gz'
    script:
        "../scripts/plots/germ_cell_cell_cycle_bar.R"

rule plot_genotype_proportions_per_germ_cluster_bar:
    input:
        sce = rules.get_adult_germ_cells.output.rds,
    output:
        pdf = 'results/plots/genotype_proportions_per_cluster_bar/genotype_proportions_per_cluster_bar.pdf',
        tsv = 'results/plots/genotype_proportions_per_cluster_bar/genotype_proportions_per_cluster_bar.tsv.gz'
    script:
        "../scripts/plots/genotype_proportions_per_germ_cluster_bar.R"

rule plot_genotype_proportions_per_celltype_bar:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
    output:
        pdf = 'results/plots/genotype_proportions_per_celltype_bar/genotype_proportions_per_celltype_bar.pdf',
        tsv = 'results/plots/genotype_proportions_per_celltype_bar/genotype_proportions_per_celltype_bar.tsv.gz'
    script:
        "../scripts/plots/genotype_proportions_per_celltype_bar.R"

rule plot_celltype_mut_vs_wt_volcano_and_ma:
    input:
        tsv = rules.adult_cmo_de.output.tsv,
    output:
        volcano = 'results/plots/celltype_mut_vs_wt_volcanos_and_ma/celltype_mut_vs_wt_volcano.pdf',
        ma = 'results/plots/celltype_mut_vs_wt_volcanos_and_ma/celltype_mut_vs_wt_ma.pdf',
        tsv = 'results/plots/celltype_mut_vs_wt_volcanos_and_ma/celltype_mut_vs_wt_volcano.tsv.gz'
    script:
        "../scripts/plots/celltype_mut_vs_wt_volcano_and_ma.R"


rule plot_germ_cell_mut_vs_wt_volcano_and_ma:
    input:
        tsv = rules.adult_cmo_germ_cell_de.output.tsv,
    output:
        volcano = 'results/plots/germ_cell_mut_vs_wt_volcanos_and_ma/germ_cell_mut_vs_wt_volcano.pdf',
        ma = 'results/plots/germ_cell_mut_vs_wt_volcanos_and_ma/germ_cell_mut_vs_wt_ma.pdf',
        tsv = 'results/plots/germ_cell_mut_vs_wt_volcanos_and_ma/germ_cell_mut_vs_wt_volcano.tsv.gz'
    script:
        "../scripts/plots/germ_cell_mut_vs_wt_volcano_and_ma.R"


rule plot_per_cell_per_sample_de_jitter_broad:
    input:
        sce = rules.find_celltypes_adult_cmo.output.rds,
        de = rules.adult_cmo_de.output.tsv,
        classifications = config.get("dfam_te_classifications"),
    output:
        pdf = 'results/plots/per_cell_per_sample_jitter/per_cell_per_sample_de_jitter.broad.pdf',
        tsv = 'results/plots/per_cell_per_sample_jitter/per_cell_per_sample_de_jitter.broad.tsv.gz'
    params:
        de_by = 'celltype',
    script:
        "../scripts/plots/per_cell_per_sample_de_jitter.R"

rule plot_per_cell_per_sample_de_jitter_germ_cells:
    input:
        sce = rules.reprocess_adult_germ_cells.output.rds,
        de = rules.adult_cmo_germ_cell_de.output.tsv,
        classifications = config.get("dfam_te_classifications"),
    output:
        pdf = 'results/plots/per_cell_per_sample_jitter/per_cell_per_sample_de_jitter.germ_cells.pdf',
        tsv = 'results/plots/per_cell_per_sample_jitter/per_cell_per_sample_de_jitter.germ_cells.tsv.gz'
    params:
        de_by = 'label',
    script:
        "../scripts/plots/per_cell_per_sample_de_jitter.R"

rule plot_de_te_heatmap_broad:
    input:
        de = rules.adult_cmo_de.output.tsv,
        sce = rules.find_celltypes_adult_cmo.output.rds,
        classifications = config.get("dfam_te_classifications"),
    output:
        pdf_fdr = 'results/plots/de_te_heatmap/de_te_heatmap.broad_fdr.pdf',
        pdf_foldchange = 'results/plots/de_te_heatmap/de_te_heatmap.broad_foldchange.pdf',
        tsv = 'results/plots/de_te_heatmap/de_te_heatmap.broad_foldchange.tsv.gz'
    params:
        de_by = 'celltype',
    script:
        "../scripts/plots/de_te_heatmap_broad.R"