# Create DGEList object of edgeR

createDGEList <- function(read_counts, group_names=rep(1,ncol(read_counts)), normalize=TRUE, norm_base=1e+6, remove_small=TRUE, verbose=FALSE, func_dir=NULL)
{
    # read_counts: a matrix of read counts and sequence length of genes at different sample conditions.
    # group_names <- c("Wildtype", "Wildtype", "Knockout", "Knockout", "Mixed", "Mixed")
    # normalize <- TRUE
    # remove_small <- TRUE
    # verbose <- FALSE

    # Load required library
    require(edgeR)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "cpb.R"), local=TRUE)
    source(file.path(func_dir, "cpb.DGEList.R"), local=TRUE)
    source(file.path(func_dir, "cpb.default.R"), local=TRUE)

    # Construct DGEList from raw read counts and comparing group names.
    dge_profile <- DGEList(counts=read_counts, group=group_names)

    # Remove very samll read counts.
    if(remove_small)
    {
        # Keep those genes with at least 1 cpb in at least the minimal number of samples in the same condition.
        kept_gene_index <- rowSums(cpb(x=dge_profile,norm.base=norm_base)>1) >= min(summary(factor(group_names)))
        dge_profile <- dge_profile[kept_gene_index,]
        # Re-compute the library sizes.
        dge_profile$samples$lib.size <- colSums(dge_profile$counts)
    }

    # Normalize read counts.
    if(normalize)
    {
        # Calculate normalization factors to scale the sizes of raw read counts.
        dge_profile <- calcNormFactors(dge_profile)
        # Print normalization factors for each sample.
        if(verbose) print(dge_profile$samples)
    }

    # Return DGEList
    return(dge_profile)
}
