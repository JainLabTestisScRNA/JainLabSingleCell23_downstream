rule import_cellranger_count:
    """
    This is a generic rule for importing some number of cellranger matrices and yielding an sce object.
    This approach aids in potentially trying different quant techniques at some point.
    it assumes raw h5 cellranger matrices (from the count workflow) and filtered matricesare available.
    """
    input: 
        raw = lambda wc: config.get("cellranger_matrices").get(wc.experiment),
        filt = lambda wc: config.get("cellranger_filtered_matrices").get(wc.experiment),
        genes = config.get("genes"),
    output: 
        rds = "results/single-cell-preproc/count/import/{experiment}.cellranger.sce.rds",
        raw = "results/single-cell-preproc/count/raw/{experiment}.cellranger.sce.rds",
    script: 
        "../scripts/single-cell-preproc/import-cellranger-count.R"

rule import_cellranger_multi:
    """
    This is a generic rule for importing some number of cellranger matrices and yielding an sce object.
    from cmo experiments. raw=all cells, rds = passing cells (per cellranger).

    In theory 'sample_assignments' could be a function of the experiment name,
    thus allowing inclusion of cmo assignments that are an output from a rule. However, this is not currently implemented.

    a raw sce is also saved to temp for use in other rules, such as ambient rna removal.
    """
    input: 
        mat = lambda wc: config.get("cellranger_matrices").get(wc.experiment),
        sample_assignments = lambda wc: config.get("cellranger_cmo_assignment").get(wc.experiment),
        genes = config.get("genes"),
    output: 
        rds = temp("results/single-cell-preproc/multi/import/{experiment}.cellranger.sce.rds"),
        raw = temp("results/single-cell-preproc/multi/raw/{experiment}.cellranger.sce.rds")
    script: 
        "../scripts/single-cell-preproc/import-cellranger.R"

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

rule preprocess_adult:
    input:
        sce = rules.import_cellranger_count.output.rds,
        raw = rules.import_cellranger_count.output.raw,
    output:
        rds = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.sce.passing_cells.rds",
        dec = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.dec.rds",
        g_dcx = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_dcx.rds", #decontx umap
        g_mito_cuts = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_mito_cuts.rds", # viz of mito cutoffs
        g_pc_elbow = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_pc_elbow.rds", # N pc choice
        g_kneeplot = "results/single-cell-preproc/preprocessed/{experiment}.cellranger.g_kneeplot.rds", # knee plots
    script:
        "../scripts/single-cell-preproc/process-adult.R"

rule integrate_adult:
    input:
        sces = expand("results/single-cell-preproc/preprocessed/{g}.cellranger.sce.passing_cells.rds", g=config.get("groupings").get("adult")),
        decs = expand("results/single-cell-preproc/preprocessed/{g}.cellranger.dec.rds", g=config.get("groupings").get("adult")),
    output:
        rds = "results/single-cell-preproc/integrated/adult.cellranger.sce.passing_cells.rds",
        g_pca_uncorrected = "results/single-cell-preproc/integrated/adult.cellranger.g_pca_uncorrected.rds",
        g_tsne_uncorrected = "results/single-cell-preproc/integrated/adult.cellranger.g_tsne_uncorrected.rds",
        g_umap_uncorrected = "results/single-cell-preproc/integrated/adult.cellranger.g_umap_uncorrected.rds",
    script:
        "../scripts/single-cell-preproc/integrate-adult.R"
