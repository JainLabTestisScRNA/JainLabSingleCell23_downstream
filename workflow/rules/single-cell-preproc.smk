rule preprocess_adult_cmo:
    """
    preproc rule specific to the adult cmo data
    """
    input:
        sce = "results/cellranger-import/{experiment}.sce.assigned.rds",
        sample_assignments = lambda wc: config.get("cellranger_cmo_assignment").get(wc.experiment),
    wildcard_constraints:
        experiment = "|".join(x for x in config.get("cellranger_cmo_assignment").keys() if 'adult' in x),
    threads: 
        4
    output:
        rds = "results/preprocessing/{experiment}.sce.passing_cells.rds",
        g_mito_cuts = "results/preprocessing/{experiment}.g_mito_cuts.rds", # viz of mito cutoffs
        g_pc_elbow = "results/preprocessing/{experiment}.g_pc_elbow.rds", # N pc choice
    script:
        "../scripts/single-cell-preproc/process-adult-cmo.R"


rule integrate_adult_cmo:
    input:
        sces = expand("results/preprocessing/{e}.sce.passing_cells.rds",e=["adult_9646_1","adult_9682_1","adult_9646_2","adult_9682_2"]),
    output:
        rds = "results/integration/adult.sce.integrated.rds",
        uncorrected = "results/integration/adult.sce.uncorrected.rds",
        mnn_out = "results/integration/adult.sce.mnn_out.rds",
    script:
        "../scripts/single-cell-preproc/integrate-adult-cmo.R"

rule cluster_adult_cmo:
    input:
        sce = "results/integration/adult.sce.integrated.rds",
    output:
        rds = "results/integration/adult.sce.integrated.clustered.rds",
        g_silhouette = "results/integration/adult.g_silhouette.rds",
        #g_coarse_clustering_sweep = "results/integration/adult.g_coarse_clustering_sweep.rds",
    script:
        "../scripts/single-cell-preproc/cluster-adult-cmo.R"

rule find_celltypes_adult_cmo:
    input:
        sce = "results/integration/adult.sce.integrated.clustered.rds",
        garnett_fl = config.get("garnett_classifier").get("broad"),
    output:
        rds = "results/integration/adult.sce.integrated.clustered.celltypes.rds",
        garnett_classifier = "results/integration/adult.garnett_classifier.rds",
        marker_check = "results/integration/adult.marker_check.rds",
    threads: 
        4
    script:
        "../scripts/single-cell-preproc/find-celltypes-adult-cmo.R"