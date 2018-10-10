# Plot heat map of gene expression distance using count per million (CPM) reads.
# The CPM trnsformation of raw read counts normalizes gene expression level
# over entire sample, so that the difference of expression level due to varied
# reading depth of each sample is eliminated and therefore different samples
# can be compared directly in terms of their gene expression profile.
# NOTE: CPM normalization doesn not change the relative expression level of
# a gene compared with that of other genes in the same sample, and therefore
# the difference of expression levels between different genes is maintained.

plotHeatMapOfGeneExpDistance <- function(gene_exp, gene_list=NULL, group_names=NULL, use_group_names=FALSE, sort_group_names=FALSE, dendrogram=FALSE, title_text=NULL, margin_bottom=5, margin_right=margin_bottom, dist_method="correlation", normalize=FALSE, log=FALSE, remove_small=FALSE, lower_limit=NULL, upper_limit=NULL, color_map="black_red_yellow_white", plot_device=NULL, plot_file=NULL, func_dir=NULL)
{
    # gene_exp: either a numerical matrix of raw read counts or a DGEList object of edgeR.
    # gene_list <- NULL # A list of selected genes to plot
    # group_names <- c("Wildtype", "Wildtype", "Knockout", "Knockout", "Mixed", "Mixed")
    # use_group_names <- FALSE
    # sort_group_names <- FALSE
    # dendrogram <- fALSE
    # title_text <- "HeatMap of Multiple Gene Samples"
    # margin_bottom <- 5
    # margin_right <- margin_bottom
    # dist_method <- "euclidean"
    # normalize <- FALSE
    # log <- FALSE
    # remove_small <- FALSE
    # lower_limit <- NULL
    # upper_limit <- NULL
    # color_map <- "black_red_yellow_white" # "blue_yellow_red" or "rainbow"
    # plot_device <- NULL # "pdf", "png", or NULL.
    # plot_file: the name of output plot file.
    # func_dir <- "/Users/granville/Documents/works/academia/MSSM/LINCS/Programs/"

    # Load required library
    require(stats)
    require(gplots)
    require(edgeR)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "createDGEList.R"), local=TRUE)
    source(file.path(func_dir, "generateColorGradient.R"), local=TRUE)
    source(file.path(func_dir, "generateRainbowHeatColors.R"), local=TRUE)

    # Determine the data type of input gene_exp.
    dge_profile <- NULL
    if(inherits(gene_exp, "DGEList"))
    {
        dge_profile <- gene_exp
        if(is.null(group_names)) group_names <- as.character(dge_profile$samples$group)
        else dge_profile$samples$group <- as.factor(group_names)
    }
    else if(is.matrix(gene_exp) && is.numeric(gene_exp))
    {
        read_counts <- gene_exp
        if(ncol(read_counts) > 1)
        {
            # Create DGEList.
            if(is.null(group_names)) group_names <- rep(1, ncol(read_counts))
            dge_profile <- createDGEList(read_counts=read_counts, group_names=group_names, normalize=normalize, remove_small=remove_small, verbose=FALSE, func_dir=func_dir)
        }
        else
        {
            # Non-fatal error.
            warning("Read counts must contain more than one observation!")
            return(NULL)
        }
    }
    else
    {
        warning("Input gene_exp must be either a numeric matrix or a DGEList object!")
        return(NULL)
    }

    # Check if a DGEList object has been created successfully.
    if(is.null(dge_profile))
    {
        warning("Failed to create DGEList object for plotting heat map!")
        return(NULL)
    }

    # Plot heat map when dge_profile is available.
    # Generate logarithmic CPM transformation of raw read counts.
    read_counts_cpm <- cpm(dge_profile, normalized.lib.sizes=normalize, log=log)
    # Extract only list genes.
    if(!is.null(gene_list))
    {
        if(is.vector(gene_list))
        {
            # Take the intersect of gene_list and read_counts_cpm due to the removal of small counts.
            gene_list_common <- intersect(gene_list, rownames(read_counts_cpm))
            read_counts_cpm <- read_counts_cpm[gene_list_common,]
        }
        else stop("Gene list must be a vector!")
    }

    # Use sample conditions to re-order the read counts matrix.
    if(use_group_names)
    {
        if(sort_group_names)
        {
            # Re-order sample columns according to their conditions.
            ss_names <- rbind(colnames(read_counts_cpm), group_names)
            ss_names <- ss_names[,sort.list(ss_names[2,], decreasing=FALSE)]
            read_counts_cpm <- read_counts_cpm[,ss_names[1,]]
            # Set column names to condition names for heat map plot.
            colnames(read_counts_cpm) <- ss_names[2,]
        }
        else
        {
            # Set column names to condition names for heat map plot.
            colnames(read_counts_cpm) <- group_names
        }
    }

    # Create a color pallette
    ncolor <- 500
    if(color_map == "black_red_yellow_white") hmcol <- generateColorGradient(pal=c("black","red","yellow","white"), n=ncolor)
    else if(color_map == "blue_yellow_red") hmcol <- generateColorGradient(pal=c("blue","yellow","red"), n=ncolor)
    else if(color_map == "rainbow") hmcol <- generateRainbowHeatColors(ncolor)
    else
    {
        warning("color_map can only be one of \"black_red_yellow_white\", \"blue_yellow_red\", or \"rainbow\"")
        return(NULL)
    }

    # Plot a heat map of DEG expression distance between multiple conditions.
    dists <- NULL
    if(dist_method != "correlation") dists <- dist(t(read_counts_cpm), method=dist_method)
    else dists <- as.dist((1-cor(read_counts_cpm))/2)
    dists_mat <- as.matrix(dists)
    rownames(dists_mat) <- colnames(dists_mat) <- colnames(read_counts_cpm)
    if(dist_method != "correlation")
    {
        if(log) key.xlab <- "logCPM"
        else key.xlab <- "CPM"
    }
    else key.xlab <- "-R*"

    # Customize color map according to lower and upper limits.
    stopifnot((is.null(lower_limit)|is.null(upper_limit)) == (is.null(lower_limit)&is.null(upper_limit)))
    breaks <- NULL
    if(!is.null(lower_limit) | !is.null(upper_limit))
    {
        dists_mat[dists_mat < lower_limit] <- NA
        dists_mat[dists_mat > upper_limit] <- NA
        breaks <- seq(lower_limit, upper_limit, length.out=(ncolor+1))
    }

    # Set cluster dendrogram options.
    rowv <- dendrogram
    symm <- TRUE
    dendrogram <- (if(dendrogram) "both" else "none")

    # Choose the plot device.
    if(!is.null(plot_device))
    {
        if(plot_device == "pdf")
        {
            plot_file_name <- paste(plot_file, plot_device, sep=".")
            pdf(file=plot_file_name, width=8.5, height=8.5, pointsize=8, onefile=TRUE)
        }
        else if(plot_device == "png")
        {
            plot_file_name <- paste(plot_file, "%02d", plot_device, sep=".")
            png(filename=plot_file_name, width=1600, height=1600, pointsize=24)
        }
        else
        {
            warning("plot_device must be one of \"pdf\", \"png\", or NULL only!")
            warning("Now set plot output to screen.")
            plot_device <- NULL
        }
    }
    else
    {
        print("Set plot output to screen.")
    }

    # Plot the heat map
    if(is.null(breaks)) heatmap.2(dists_mat, Rowv=rowv, symm=TRUE, dendrogram=dendrogram, scale="none", trace="none", col=hmcol, key.title="Expression Histogram", key.xlab=key.xlab, keysize=1, margin=c(margin_bottom,margin_right), main=title_text)
    else heatmap.2(dists_mat, Rowv=rowv, symm=TRUE, dendrogram=dendrogram, scale="none", trace="none", col=hmcol, key.title="Expression Histogram", key.xlab=key.xlab, breaks=breaks, keysize=1, margin=c(margin_bottom,margin_right), main=title_text)

    # Close the PDF plot.
    if(!is.null(plot_device)) dev.off()

    # Return CPM read counts implicitly.
    invisible(read_counts_cpm)
}
