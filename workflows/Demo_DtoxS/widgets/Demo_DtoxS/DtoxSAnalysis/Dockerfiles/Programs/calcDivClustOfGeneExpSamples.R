# Calculate divisive hierarchical clusters of given gene expression samples.

calcDivClustOfGeneExpSamples <- function(read_counts, metric="correlation")
{
    # read_counts: a matrix of gene expression samples.
    # metric: "correlation", "euclidean", or "manhattan".
    
    # Load needed functions.
    require(stats)
    require(cluster)
    
    # Check argument validity
    if(!is.matrix(read_counts) || !is.numeric(read_counts))
    {
        warning("Input argument must be a numeric matrix!")
        return(NULL)
    }
    
    # Check the number of samples.
    if(ncol(read_counts) < 2)
    {
        warning("At least 2 samples are required for calculate divisive clusters!")
        return(NULL)
    }
    
    # Decide which method to use for calculation.
    if(metric == "correlation")
    {
        # Calculate the distance between samples based on their correlation.
        dist_mat <- as.dist((1-cor(read_counts))/2)
        # Calculate divisive clusters.
        div_clust <- diana(dist_mat, metric=metric)
    }
    else div_clust <- diana(t(read_counts))
    
    # Return the results.
    return(div_clust)
}
