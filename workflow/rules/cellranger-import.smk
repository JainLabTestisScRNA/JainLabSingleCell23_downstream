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
        rds = "results/cellranger-import/{experiment}.sce.assigned.rds",
        raw = "results/cellranger-import/{experiment}.sce.raw.rds"
    script: 
        "../scripts/single-cell-preproc/import-cellranger.R"
