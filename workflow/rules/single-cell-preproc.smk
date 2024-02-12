"""
Spec for preprocessed sce objects.
These are the objects that will be used for all downstream analysis.

sce with coldata including the following columns:
label - coarse cluster
sublabel - fine cluster
celltype  - celltype at the sublabel level, low uniformity is allowed (e.g. ie some major clusters are mixed from sublabels of slighlty different celltypes)
macro_celltype - celltype at the macro level (e.g. ie pretty uniform in major clusters)
label2 - coarse cluster + macro_celltype
sublabel2 - fine cluster + macro_celltype
Sample - genotype (e.g. cmo if that is the data in question)
"""

rule preprocess_juvenile:
    """
    preproc rule specific to the juvenile cmo data
    """
    input:
        sce = "results/single-cell-preproc/multi/import/juvenile_13d_wt_null.cellranger.sce.rds",
        raw = "results/single-cell-preproc/multi/raw/juvenile_13d_wt_null.cellranger.sce.rds",
        celltype_calls = lambda wc: config.get("celltype_assignment").get("juvenile_13d_wt_null"),
        sample_assignments = lambda wc: config.get("cellranger_cmo_assignment").get("juvenile_13d_wt_null"),
    output:
        rds = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.passing_cells.rds",
        g_dcx = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.g_dcx.rds", #decontx umap
        g_mito_cuts = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.g_mito_cuts.rds", # viz of mito cutoffs
        g_pc_elbow = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.g_pc_elbow.rds", # N pc choice
        g_coarse_clustering_sweep = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.g_coarse_clustering_sweep.rds", # viz of clustering sweeps
        nn_clust = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.nn_clust.rds", # igraph obj of clustering
        g_silhouette = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.g_silhouette.rds", # silhouette plot
        sce_mito_warning = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.mito_warning.rds", # cells below hard cutoff, but still very high mito
        sce_sertoli = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.sertoli.rds", # sertoli cells
        sce_germ_cell = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.germ_cell.rds", # germ cells
        sce_somatic_non_macro = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.somatic_non_macro.rds", # somatic cells that are not macrophages
        sce_macrophage = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.macrophage.rds", # macrophages
    script:
        "../scripts/single-cell-preproc/process-juvenile.R"

rule preprocess_adult_cmo:
    """
    preproc rule specific to the adult cmo data
    """
    input:
        sce = "results/single-cell-preproc/multi/import/{experiment}.cellranger.sce.rds",
        raw = "results/single-cell-preproc/multi/raw/{experiment}.cellranger.sce.rds",
        sample_assignments = lambda wc: config.get("cellranger_cmo_assignment").get(wc.experiment),
        garnett_fl = lambda wc: config.get("garnett_classifier").get(wc.experiment),
    wildcard_constraints:
        experiment = "|".join(x for x in config.get("cellranger_cmo_assignment").keys() if 'adult' in x),
    threads: 
        8
    output:
        rds = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.sce.passing_cells.rds",
        g_dcx = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_dcx.rds", #decontx umap
        g_mito_cuts = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_mito_cuts.rds", # viz of mito cutoffs
        g_pc_elbow = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_pc_elbow.rds", # N pc choice
        g_coarse_clustering_sweep = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_coarse_clustering_sweep.rds", # viz of clustering sweeps
        nn_clust = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.nn_clust.rds", # igraph obj of clustering
        g_silhouette = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_silhouette.rds", # silhouette plot
        sce_mito_warning = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.sce.mito_warning.rds", # cells below hard cutoff, but still very high mito
        sce_sertoli = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.sce.sertoli.rds", # sertoli cells
        sce_germ_cell = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.sce.germ_cell.rds", # germ cells
        sce_somatic_non_macro = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.sce.somatic_non_macro.rds", # somatic cells that are not macrophages
        sce_macrophage = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.sce.macrophage.rds", # macrophages
    script:
        "../scripts/single-cell-preproc/process-adult-cmo.R"


rule integrate_adult_cmo:
    input:
        sces = expand("results/single-cell-preproc/preprocessed/{e}.cellranger.sce.passing_cells.rds",e=["adult_9646_combined","adult_9682_combined"]),
        garnett_fl = config.get("garnett_classifier").get("adult"),
    output:
        rds = "results/single-cell-preproc/integrated/adult.cellranger.sce.integrated.rds",
        uncorrected = "results/single-cell-preproc/integrated/adult.cellranger.sce.uncorrected.rds",
        garnett_classifier = "results/single-cell-preproc/integrated/adult.cellranger.garnett_classifier.rds",
        g_silhouette = "results/single-cell-preproc/integrated/adult.cellranger.g_silhouette.rds",
        g_coarse_clustering_sweep = "results/single-cell-preproc/integrated/adult.cellranger.g_coarse_clustering_sweep.rds",
        marker_check = "results/single-cell-preproc/integrated/adult.cellranger.marker_check.rds",
        mnn_out = "results/single-cell-preproc/integrated/adult.cellranger.sce.mnn_out.rds",
    script:
        "../scripts/single-cell-preproc/integrate-adult-cmo.R"