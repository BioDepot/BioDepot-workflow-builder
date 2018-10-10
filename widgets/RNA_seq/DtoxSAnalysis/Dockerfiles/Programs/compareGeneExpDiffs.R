# Compare the pairwise differences between multiple gene expression profiles.

compareGeneExpDiffs <- function(read_counts, group_names, group_pairs, dispersion="auto", min_rep=3, rep_req="all", fdr=0.05, deg_only=TRUE, normalize=TRUE, norm_base=1e+6, remove_small=TRUE, exp_info=NULL, deg_read_counts_file=NULL, plot_bcv=FALSE, plot_smear=FALSE, deg_plots_file=NULL, pt_chars=c(16,19), sub_char=FALSE, r_compat=TRUE, verbose=FALSE, func_dir=NULL)
{
    # read_counts: a matrix of read counts and sequence length of genes at different sample conditions.
    # group_names <- c("Wildtype", "Wildtype", "Knockout"", "Knockout"", "Mixed", "Mixed")
    # group_pairs <- c("Wildtype,Knockout", "Wildtype,Mixed", "Knockout,Mixed") (group_pairs can be NULL) # Pairs of the State fields in targets.
    # dispersion <- "auto"
    # min_rep <- 3 # The minimum number of sample replicates of the groups to be compared.
    # rep_req <- "all" or "any" # A minimum replicate number for all groups or for at least one of them.
    # fdr <- 0.05
    # deg_only <- TRUE
    # normalize <- TRUE
    # remove_small <- TRUE
    # exp_info <- NULL # Experimental information including drug group, cell line, time point, etc.
    # deg_read_counts_file <- "Human.A-Hour.48-Plate.1-Calc-CTRL.TRS+LOP.tsv"
    # plot_bcv <- FALSE # Don't plot BCV.
    # plot_smear <- FALSE # Don't plot mean-variance.
    # deg_plots_file <- "Human-48-CTRL-TRS-DEG-Plots-Exp3.all.refseq.umi.0.1.pdf"
    # pt_chars <- c(16,19) # Two kinds of point characters for plotBCV and plotSmear.
    # sub_char <- FALSE # Flag of replacing the unsuitable characters in state names with underscore, so that the state names can be used as variable names.
    # r_compat <- TRUE # Generate R-compatible or view-compatible data file.
    # verbose <- FALSE
    # func_dir <- "/Users/granville/Documents/works/academia/MSSM/LINCS/Programs/"

    # Load required library
    require(tools)
    require(stats)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "compareGeneExpDiff.R"), local=TRUE)

    # Initialize the list of DEG read counts.
    deg_read_counts_list <- list()

    # Find corresponding DEGs in each comparison pair.
    # Traverse through all pairs for comparison.
    for(pair_idx in 1:nrow(group_pairs))
    {
        # Get a pair of group names for a list of sample replicates to compare.
        group_pair <- as.character(group_pairs[pair_idx,])

        # Print comparison information to console.
        if(verbose) print(paste0(group_pair[1], " vs ", group_pair[2], ":"))

        # Prepare the dataset of current sample pair for comparison.
        group_pair_idx <- group_names==as.character(group_pair)[1]|group_names==as.character(group_pair)[2]
        group_names_compare <- group_names[group_pair_idx]
        read_counts_compare <- read_counts[,group_pair_idx]
        # Remove the rows containing any NA.
        read_counts_compare <- read_counts_compare[rowSums(is.na(read_counts_compare))==0,]
        # Remove the rows containing all zero read counts.
        read_counts_compare <- read_counts_compare[rowSums(read_counts_compare)>0,]
        sample_names_compare <- colnames(read_counts_compare)

        # Compare the expression profile of current pair of samples.
        deg_counts_stats <- compareGeneExpDiff(read_counts=read_counts_compare, group_names=group_names_compare, group_pair=group_pair, dispersion=dispersion, min_rep=min_rep, rep_req=rep_req, fdr=fdr, deg_only=deg_only, normalize=normalize, norm_base=norm_base, remove_small=remove_small, exp_info=exp_info, deg_read_counts_file=deg_read_counts_file, plot_bcv=plot_bcv, plot_smear=plot_smear, deg_plots_file=deg_plots_file, pt_chars=pt_chars, r_compat=r_compat, verbose=verbose, func_dir=func_dir)

        # Replcae "+" character to "_" for using as drug names as the names of list elements.
        group_pair_1 <- group_pair[1]
        group_pair_2 <- group_pair[2]
        if(sub_char)
        {
            group_pair_1 <- gsub("\\+", "_", group_pair_1)
            group_pair_2 <- gsub("\\+", "_", group_pair_2)
        }
        # Save current DEG read counts and statistics to the list
        deg_read_counts_list[[paste(group_pair_1,group_pair_2,sep=".")]] <- deg_counts_stats
    }

    # Return the list of DEG read counts.
    return(deg_read_counts_list)
}
