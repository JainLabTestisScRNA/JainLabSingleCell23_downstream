"""
As much as I hate to have specific rules for each trajectory, there is a lot of tuning involved in each of them,
so will have a separate rule for each trajectory.
"""
rule juvenile_germ_trajectory:
    input:
        sce = "results/single-cell-preproc/preprocessed/juvenile_13d_wt_null.cellranger.sce.germ_cell.rds",
    output:
        rds = "results/trajectory/juvenile_13d_wt_null.cellranger.sce.germ_cell.trajectory.rds",
    script:
        "../scripts/trajectory/juvenile_germ_trajectory.R"

rule adult_germ_trajectory:
    input:
        sce = "results/single-cell-preproc/integrated/adult.cellranger.sce.integrated.rds",
    output:
        rds = "results/trajectory/adult.cellranger.sce.germ_cell.trajectory.rds",
    script:
        "../scripts/trajectory/adult_germ_trajectory.R"



"""
Thankfully I think this requires low/no tuning, so these can be general rules
"""
rule within_experiment_subset_tradeseq:
    """
    Assumes slingshot has been run on the input.
    """
    input:
        sce = "results/trajectory/{experiment}.{quant_tool}.sce.{subset}.trajectory.rds",
    threads: 2,
    params:
        knots = 5,
    output:
        rds = "results/tradeseq/{experiment}.{quant_tool}.sce.{subset}.tradeseq.rds",
    script:
        "../scripts/trajectory/within_experiment_subset_tradeseq.R"

rule inter_condition_deg_tradeseq:
    """
    Assumes tradeseq has been run on the input.
    """
    input:
        sce = rules.within_experiment_subset_tradeseq.output.rds,
    output:
        deg_rds = "results/tradeseq/{experiment}.{quant_tool}.{subset}.deg.rds",
    script:
        "../scripts/trajectory/inter_condition_deg_tradeseq.R"

rule inter_condition_gene_clusters_tradeseq:
    """
    Assumes tradeseq has been run on the input. Runs only on signficant DEGs.
    """
    input:
        sce = rules.within_experiment_subset_tradeseq.output.rds,
        deg = rules.inter_condition_deg_tradeseq.output.deg_rds,
    output:
        gene_clusters_tsv = "results/tradeseq/{experiment}.{quant_tool}.{subset}.gene_clusters.tsv",
        smoothed_cluster_timecourse_rds = "results/tradeseq/{experiment}.{quant_tool}.{subset}.smoothed_cluster_timecourse.rds",
        g_k_choice = "results/tradeseq/{experiment}.{quant_tool}.{subset}.g_k_choice.rds",
    script:
        "../scripts/trajectory/inter_condition_gene_clusters_tradeseq.R"

rule inter_condition_gene_clusters_clusterprofiler:
    """
    Runs only on signficant DEGs. Takes input from tradeseq->conditionTest->clustering workflow.
    uses expressed genes as 'universe' - see script for full criteria for 'expressed'
    """
    input:
        gene_clusters_tsv = rules.inter_condition_gene_clusters_tradeseq.output.gene_clusters_tsv,
        smoothed = rules.inter_condition_gene_clusters_tradeseq.output.smoothed_cluster_timecourse_rds,
        sce = rules.within_experiment_subset_tradeseq.output.rds,
    output:
        #deg_enrichment_tsv = "results/tradeseq/{experiment}.{quant_tool}.{subset}.deg_cluster_enrichment.tsv",
        clusterprofiler_rds = "results/tradeseq/{experiment}.{quant_tool}.{subset}.clusterprofiler.rds",
    script:
        "../scripts/trajectory/inter_condition_gene_clusters_clusterprofiler.R"

rule dysregulated_te_seqs:
    """
    Input is table of DE features from tradeseq. Output is fasta of sequences of dysregulated TEs.
    """
    input:
        deg = rules.inter_condition_deg_tradeseq.output.deg_rds,
        fasta = config.get("te_sequences")
    output:
        fa = "results/tradeseq/{experiment}.{quant_tool}.{subset}.dysregulated_tes.fasta",
    script:
        "../scripts/trajectory/dysregulated_te_seqs.R"   