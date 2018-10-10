# Calculate and plot divisive clusters of gene expression read counts.

clusterGeneExpSamples <- function(read_counts, plot_clust=FALSE, drug_name=NULL, cell_line=NULL, plate_number=NULL, dist_cutoffs=NULL, dist_cutoff_outlier=0.01, dist_cutoff_group=0.015, min_samples=3, ylim=NULL, hline=NULL, title_text=NULL, func_dir=NULL)
{
    # Load required library
    require(stats)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "getCutoff.R"), local=TRUE)
    source(file.path(func_dir, "calcInvExp10ofDecimal.R"), local=TRUE)
    source(file.path(func_dir, "calcDivClustOfGeneExpSamples.R"), local=TRUE)
    source(file.path(func_dir, "plotDivClustOfGeneExpSamples.R"), local=TRUE)

    # Check argument validity
    if(!is.matrix(read_counts) || !is.numeric(read_counts))
    {
        warning("Input argument must be a numeric matrix!")
        return(NULL)
    }

    # Check the availability of the arguments related to cluster plot.
    if(plot_clust)
    {
        if(is.null(drug_name) || is.null(cell_line) || is.null(plate_number) || is.null(dist_cutoffs))
        {
            warning("Input arguments (drug_name, cell_line, plate_number, and dist_cutoffs) must be set!")
            return(NULL)
        }
    }

    # Determine the availability of input arguments: ylim and hline.
    if(is.null(ylim)) return_ylim <- TRUE
    else return_ylim <- FALSE
    if(is.null(hline)) return_hline <- TRUE
    else return_hline <- FALSE

    # Calculate divisive clusters for the samples in current cell line and drug-treated condition.
    divclust <- calcDivClustOfGeneExpSamples(read_counts=read_counts, metric="correlation")

    # Plot divisive clusters of orginal samples.
    if(plot_clust && is.matrix(read_counts))
    {
        # Calculate the upper limit of Y axis for the cluster tree plot.
        if(is.null(ylim))
        {
            dist_max <- (1 - min(cor(read_counts))) / 2
            inv_exp10 <- calcInvExp10ofDecimal(dist_max)
            dist_max_mag <- dist_max * 10^inv_exp10
            dist_max_mag_round <- round(dist_max_mag)
            dist_max_mag_floor <- floor(dist_max_mag)
            dist_max_mag_ceiling <- ceiling(dist_max_mag)
            if(dist_max_mag_round == dist_max_mag_floor) ylim <- (dist_max_mag_floor + 0.5) / 10^inv_exp10
            else ylim <- dist_max_mag_ceiling / 10^inv_exp10
        }
        # Calculate the cutoff line for the cluster tree plot.
        if(is.null(hline))
        {
            # Set cutoff values for outlier samples.
            cutoff <- getCutoff(state=drug_name, cell=cell_line, plate=plate_number, cutoffs=dist_cutoffs, single=dist_cutoff_outlier, group=dist_cutoff_group)
            if(ncol(read_counts) > min_samples) hline <- cutoff[1]
            else hline <- cutoff[2]
        }
        # Adjust ylim if hline is above it.
        if(is.infinite(hline)) hline <- NULL
        if(is.null(ylim))
        {
            if(!is.null(hline)) ylim <- ceiling(hline * 10^inv_exp10) / 10^inv_exp10
        }
        else
        {
            if(!is.null(hline) && hline > ylim) ylim <- ceiling(hline * 10^inv_exp10) / 10^inv_exp10
        }
        if(!is.null(ylim)) ylim <- c(0, ylim)
        # Plot divisive clusters of filtered samples.
        if(is.null(title_text)) title_text <- paste(drug_name, "on Cell", cell_line, "in Plate", plate_number)
        plotDivClustOfGeneExpSamples(divclust=divclust, main=title_text, hline=hline, ylim=ylim, cex.main=1.75)
        if(!is.null(ylim)) ylim <- ylim[2]
    }

    # Return cluster results.
    if(!is.null(divclust))
    {
        if(!return_ylim && !return_hline) result <- divclust
        else if(return_ylim && !return_hline)
        {
            if(!is.null(ylim)) result <- list(clust=divclust, ylim=ylim)
            else result <- divclust
        }
        else if(!return_ylim && return_hline) result <- list(clust=divclust, hline=hline)
        else
        {
            if(!is.null(ylim)) result <- list(clust=divclust, ylim=ylim, hline=hline)
            else result <- list(clust=divclust, hline=hline)
        }
    }
    else
    {
        if(!return_ylim && !return_hline) result <- NULL
        else if(return_ylim && !return_hline)
        {
            if(!is.null(ylim)) result <- ylim
            else result <- NULL
        }
        else if(!return_ylim && return_hline) result <- hline
        else
        {
            if(!is.null(ylim)) result <- list(ylim=ylim, hline=hline)
            else result <- hline
        }
    }
    invisible(result)
}
