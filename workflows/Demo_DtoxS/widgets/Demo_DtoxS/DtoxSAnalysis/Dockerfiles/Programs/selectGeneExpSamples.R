# Select the sample replcates of gene expression data according to given cutoff values.

selectGeneExpSamples <- function(gene_exp, cutoff_outlier=0.01, cutoff_group=cutoff_outlier, min_samples=3, keep_under_samples=FALSE, verbose=FALSE, func_dir=NULL)
{
    # gene_exp
    # A data matrix of gene expression samples or a diana object.
    #
    # cutoff_outlier: 0.01
    # The distance cutoff for removing outliers of a sample group with a group
    # size greater than min_samples replicates.
    #
    # cutoff_group: cutoff_outlier
    # The distance cutoff for discarding entire sample group of a group size
    # equal to min_samples replicates, such that entire sample group will be
    # kept if its size is greater than cutoff_outlier, but less than cutoff_group.
    # This threshold is also used for a sample group with size less than
    # min_samples replicates if keep_under_samples is TRUE in order to retain
    # remaining sample replicates after cut off for later combination with
    # other condition samples.
    #
    # min_samples: 3
    # The minimum number of sample replicates for determining which cutoff to
    # use: cutoff_outlier or cutoff_group.
    #
    # keep_under_samples: FALSE
    # The flag indicating whether or not to keep the remaining samples with
    # the number of replicates is less than min_samples after cut off. Such
    # case can happen regardless of the replicate number of original samples.
    #
    # verbose: FALSE
    # The switch to turn on printing debug information.
    #
    # func_dir: "/Users/granville/Documents/Repos/works/MSSM/LINCS/Programs/Library"
    # The directory for loading user-defined functions used.

    # Load required library
    source(file.path(func_dir, "calcDivClustOfGeneExpSamples.R"), local=TRUE)
    source(file.path(func_dir, "lappend.R"), local=TRUE)

    # Check input minimum number of samples.
    if(min_samples < 2)
    {
        warning("Input min_samples must be 2 at least!")
        return(NULL)
    }

    # Check input gene expression data.
    divclust <- NULL
    if(inherits(gene_exp, "diana")) divclust <- gene_exp
    else if(inherits(gene_exp, "matrix"))
    {
        # Calculate divisive hierarchical clusters of given gene expression samples.
        divclust <- calcDivClustOfGeneExpSamples(gene_exp, metric="correlation")
    }
    else
    {
        warning("Input gene_exp must be either a numeric matrix or a diana object!")
        return(NULL)
    }

    # Initialize the variables for selection.
    sample_names_sel <- NULL
    sample_height_sel <- NULL

    # Select the sample replicates whose distance is within preset cutoffs.
    if(!is.null(divclust))
    {
        sample_number <- length(divclust$order)
        divclust_height_max <- max(divclust$height)
        if((sample_number > min_samples && divclust_height_max > cutoff_outlier) || (sample_number <= min_samples && divclust_height_max > cutoff_group))
        {
            # Use cutoff_outlier or cutoff_group thresholds to cut off a sample
            # group whose replicate number is greather than or not greater than
            # min_samples.

            # Only call the "cut" function if the maximum height of the original
            # cluster is greater than cutoff_outlier, because the "cut" function will
            # cut the original cluster into upper and lower parts no matter its maximum
            # height is greater than the given cutoff_outlier or not.

            # cutoff_outlier is used for sample groups with more than or equal to 4 replicates.
            divdend <- as.dendrogram(divclust)
            # Use cutoff_outlier for a sample group with more than min_samples replicates.
            # Use cutoff_group for a sample group with equal to or less than min_samples replicates.
            if(sample_number > min_samples) divdend_divided <- cut(divdend, cutoff_outlier)
            else divdend_divided <- cut(divdend, cutoff_group)
            branch_heights <- NULL
            leaf_numbers <- NULL
            branch_leaves <- list()
            for(divdend_branch in divdend_divided$lower)
            {
                if(verbose) print(divdend_branch)
                if(!is.leaf(divdend_branch))
                {
                    # Sub-cluster
                    branch_heights <- c(branch_heights, attr(divdend_branch, "height"))
                    leaf_numbers <- c(leaf_numbers, nobs(divdend_branch))
                    branch_leaves <- lappend(branch_leaves, labels(divdend_branch))
                    if(verbose) print(paste("-->", nobs(divdend_branch), "Leaves:", paste0(labels(divdend_branch), collapse=", ")))
                }
            }

            # Find (multiple) sub-clusters with the same most number of leaf nodes
            # among all sub-clusters.
            leaf_number_max <- max(leaf_numbers)
            leaf_number_max_idx <- leaf_numbers==leaf_number_max
            # Find the sub-clusters with the least height (cluster diameter) among
            # those sub-clusters with the same most number of leaf nodes.
            branch_heights_max_leaves <- branch_heights[leaf_number_max_idx]
            branch_height_min <- min(branch_heights_max_leaves)
            branch_height_min_idx <- branch_heights==branch_height_min
            branch_leaves_max_leaf_min_height <- unlist(branch_leaves[branch_height_min_idx])

            # Decide whether or not to select, discard, or keep the remaining cluster.
            n_leaves <- length(branch_leaves_max_leaf_min_height)
            if(n_leaves >= min_samples || (n_leaves < min_samples && n_leaves > 0 && keep_under_samples))
            {
                # A minimum of min_samples sample replicates are required.
                sample_names_sel <- branch_leaves_max_leaf_min_height
                sample_height_sel <- branch_height_min
            }
            else warning("This sample group is discarded!")
        }
        else
        {
            # When a sample group is below threshold, keep it if:
            # 1) This sample group has no less than min_samples replicates.
            # 2) This sample group has less than min_samples replicates, but
            #    keep_under_samples is TRUE.
            if((sample_number > min_samples && divclust_height_max <= cutoff_outlier) || (sample_number == min_samples && divclust_height_max <= cutoff_group) || (sample_number < min_samples && divclust_height_max <= cutoff_group && keep_under_samples))
            {
                sample_names_sel <- divclust$order.lab[sort.list(divclust$order,decreasing=FALSE)]
                sample_height_sel <- divclust_height_max
            }
            else warning("This sample group is discarded!")
        }
    }
    else
    {
        # In case when gene_exp is a matrix and calcDivClustOfGeneExpSamples returns NULL,
        # then return column names of gene_exp if keep_under_samples is TRUE.
        if(keep_under_samples) sample_names_sel <- colnames(gene_exp)
    }

    # The selected sample replicates
    if(!is.null(sample_names_sel))
    {
        if(verbose) print("The selected sample replicates:")
        if(verbose) print(sample_names_sel)
        if(verbose) print(paste("And the height is", sample_height_sel))
    }
    else
    {
        if(verbose) print("No sample replicate is selected")
    }

    # Return the selected sample replicates
    return(sample_names_sel)
}
