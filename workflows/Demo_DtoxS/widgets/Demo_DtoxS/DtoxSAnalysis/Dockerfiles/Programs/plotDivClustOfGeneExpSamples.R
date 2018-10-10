# Plot divisive clusters for the dataset given state name and cell line.

plotDivClustOfGeneExpSamples <- function(divclust, main=NULL, ylim=NULL, hline=NULL, ...)
{
    # Load needed functions.
    require(cluster)
    
    # Check the data type.
    if(!inherits(divclust, "diana"))
    {
        warning("The input argument must be an object of \"diana\" class!")
        invisible(NULL)
    }
    
    # Plot divisive clusters.
    if(length(divclust$order) > 2)
    {
        if(is.null(ylim)) plot(divclust, main=main, xlab="Sample", which.plots=2, ask=FALSE, ...)
        else plot(as.dendrogram(divclust), main=main, ylab="Height", ylim=ylim, ...)
        # Add a horizontal line.
        if(!is.null(hline)) abline(h=hline, col="deeppink", lty=4)
    }
    else warning("At least 3 samples are required for cluster plot!")
    
    # Return the results.
    invisible(divclust)
}
