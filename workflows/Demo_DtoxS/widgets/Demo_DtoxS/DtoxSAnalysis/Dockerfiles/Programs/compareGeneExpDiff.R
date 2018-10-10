# Compare the difference between two gene expression profiles.

compareGeneExpDiff <- function(read_counts, group_names, group_pair, dispersion="auto", min_rep=3, rep_req="all", fdr=0.05, deg_only=TRUE, normalize=TRUE, norm_base=1e+6, remove_small=TRUE, exp_info=NULL, deg_read_counts_file=NULL, plot_bcv=FALSE, plot_smear=FALSE, deg_plots_file=NULL, pt_chars=c(16,19), r_compat=TRUE, verbose=FALSE, func_dir=NULL)
{
    # read_counts: a matrix of read counts and sequence length of genes at different sample conditions.
    # group_names <- c("Wildtype", "Wildtype", "Knockout"", "Knockout"", "Mixed", "Mixed")
    # group_pair <- c("Wildtype","Knockout") # A pair of the State fields in targets.
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
    # deg_plots_file <- "Human.A-Hour.48-Plate.1-DEG.0.1-Plots-CTRL.TRS+LOP.pdf"
    # pt_chars <- c(16,19) # Two kinds of point characters for plotBCV and plotSmear.
    # r_compat <- TRUE # Generate R-compatible or view-compatible data file.
    # verbose <- FALSE
    # func_dir <- "/Users/granville/Documents/works/academia/MSSM/LINCS/Programs/"

    # Load required library
    require(tools)
    require(stats)

    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "format_decimal.R"), local=TRUE)
    source(file.path(func_dir, "createDGEList.R"), local=TRUE)

    # Check the arguments of read_counts and group_names.
    if(ncol(read_counts) != length(group_names))
    {
        warning("The column number of read_counts must be equal to the length of group_names!")
        return(NULL)
    }

    # Check the pair argument.
    if(length(group_pair) != 2)
    {
        warning("Group pair must be a character vector of length 2!")
        return(NULL)
    }
    if(length(unique(group_pair)) < length(group_pair))
    {
        warning("Group pair cannot contain identical elements!")
        return(NULL)
    }

    # Check the value of rep_req.
    rep_req <- tolower(rep_req)
    if(rep_req!="all" && rep_req!="any")
    {
        warning("rep_req must be either \"all\" or \"any\"!")
        return(NULL)
    }

    # Extract the main and extension names from output file of DEG read counts.
    if(!is.null(deg_read_counts_file))
    {
        deg_read_counts_file_base <- basename(deg_read_counts_file)
        deg_read_counts_file_dir <- dirname(deg_read_counts_file)
        deg_read_counts_file_main_name <- file_path_sans_ext(deg_read_counts_file_base)
        deg_read_counts_file_ext_name <- file_ext(deg_read_counts_file_base)
    }
    else
    {
        deg_read_counts_file_main_name <- NULL
        deg_read_counts_file_ext_name <- NULL
    }
    # Generate the file name for the dataset of current DEG read counts.
    if(!is.null(deg_read_counts_file_main_name))
    {
        # Customize file name.
        deg_read_counts_file_name <- paste(deg_read_counts_file_main_name, paste(group_pair[1],group_pair[2],sep="."), sep="-")
        deg_read_counts_file_name <- paste(deg_read_counts_file_name, deg_read_counts_file_ext_name, sep=".")
        # Add directory path.
        deg_read_counts_file_name <- file.path(deg_read_counts_file_dir, deg_read_counts_file_name)
    }
    else deg_read_counts_file_name <- NULL

    # Extract the main and extension names from output file of DEG plots.
    if(!is.null(deg_plots_file))
    {
        deg_plots_file_base <- basename(deg_plots_file)
        deg_plots_file_dir <- dirname(deg_plots_file)
        deg_plots_file_main_name <- file_path_sans_ext(deg_plots_file_base)
        deg_plots_file_ext_name <- file_ext(deg_plots_file_base)
    }
    else
    {
        deg_plots_file_main_name <- NULL
        deg_plots_file_ext_name <- NULL
    }
    # Generate the file name for the plots of current read counts and DEGs.
    if(!is.null(deg_plots_file_main_name))
    {
        # Customize file name.
        deg_plots_file_name <- paste(deg_plots_file_main_name, paste(group_pair[1],group_pair[2],sep="."), sep="-")
        deg_plots_file_name <- paste(deg_plots_file_name, deg_plots_file_ext_name, sep=".")
        # Add directory path.
        deg_plots_file_name <- file.path(deg_plots_file_dir, deg_plots_file_name)
    }
    else deg_plots_file_name <- NULL

    # Generate the output plot PDF file.
    if(!is.null(deg_plots_file_name) && deg_plots_file_name!="" && (plot_bcv || plot_smear))
    {
        # Open a PDF file in letter size to plot
        pdf(deg_plots_file_name, width=8.5, height=8.5, onefile=TRUE)
        plot_pdf_flag <- TRUE
    }
    else plot_pdf_flag <- FALSE

    # Only include the samples to be compared.
    compare_idx <- rep(FALSE, length(group_names))
    for(item in group_pair) compare_idx <- (compare_idx | (group_names==item))
    group_names_compare <- group_names[compare_idx]
    group_names_number <- length(unique(group_names_compare))
    if(group_names_number < length(group_pair))
    {
        warning("The samples to compare cannot be found!")
        # Close plot device
        if(plot_pdf_flag)
        {
            dev.off()
            plot_pdf_flag = FALSE
            file.remove(deg_plots_file_name)
        }
        return(NULL)
    }
    sample_names_compare <- colnames(read_counts)[compare_idx]
    read_counts_compare <- read_counts[,sample_names_compare]
    read_counts_compare <- read_counts_compare[rowSums(is.na(read_counts_compare))==0,]

    # Determine if every sample group has only one replicate.
    # any: requires replicate number for each group must be greater than the threshold.
    # all: requires replicate number for at least one group must be greater than the threshold.
    rep_req_error <- table(group_names_compare) < min_rep
    if(if(rep_req=="all") any(rep_req_error) else all(rep_req_error))
    {
        warning(paste0(paste("A minimum of", min_rep, "replicates are required for", (if(rep_req=="all") "all groups" else "at least one group")), "!"))

        # Close plot device
        if(plot_pdf_flag)
        {
            dev.off()
            plot_pdf_flag = FALSE
            file.remove(deg_plots_file_name)
        }

        return(NULL)
    }

    # Construct DGEList from raw read counts and comparing group names.
    dge_profile <- createDGEList(read_counts=read_counts_compare, group_names=group_names_compare, normalize=normalize, remove_small=remove_small, norm_base=norm_base, verbose=verbose, func_dir=func_dir)

    # Estimate common, trended, and tagwise negative Binomial dispersions for entire dataset containing all combinations of interested factors.
    dge_profile <- estimateDisp(dge_profile)
    # Calculate common BCV when there are more than one sample.
    if(verbose) print(paste("BCV =", format_decimal(sqrt(dge_profile$common.dispersion))))
    # Plot biological coefficient of variation
    if(plot_bcv)
    {
        title_text <- paste("Common BCV of", paste0(group_pair, collapse="/"))
        if(!is.null(exp_info) && !is.na(exp_info["cell"])) title_text <- paste(title_text, exp_info["cell"])
        if(!is.null(exp_info) && !is.na(exp_info["time"])) title_text <- paste(title_text, "in", exp_info["time"])
        title_text <- paste(title_text, "is", format_decimal(sqrt(dge_profile$common.dispersion)))
        plotBCV(dge_profile, xlab="Average logCPM", ylab="Biological Coefficient of Variation", main=title_text, pch=pt_chars[1], cex.main=1.75)
    }

    # Print comparison information to console.
    if(verbose) print(paste0(group_pair[1], " vs ", group_pair[2], ":"))

    # Perform exact tests for gene-wise differences between two groups
    # of negative-binomially distributed counts, and obtain multiple
    # statistical measures for each gene, e.g. logFC, logCPM, and p-value.
    dge_diffs <- exactTest(dge_profile, pair=group_pair, dispersion=dispersion)
    # Extract statistic measures and calculate adjusted p-values for each genes.
    # Do not sort top_degs, so that it can be spliced with read counts table.
    top_degs <- topTags(dge_diffs, n=nrow(dge_diffs$table), sort.by="none")
    # Print top 10 DEGs ranked by their p-value or absolute log-fold change.
    if(verbose) print((top_degs$table[sort.list(top_degs$table[,"FDR"], decreasing=FALSE),])[1:10,])

    # Merge raw read counts, normalized read counts, read counts statistics,
    # and DEG statistics of all genes into one matrix.
    read_counts_norm <- cpm(dge_profile)
    colnames(read_counts_norm) <- paste(colnames(read_counts_norm), "Norm", sep=".")
    dge_counts_stats <- cbind(dge_profile$counts, read_counts_norm, top_degs$table)

    # Classify DEG statistics as up- or down-regulation, or not significant.
    dge_pattern <- decideTestsDGE(dge_diffs, p.value=fdr)
    if(verbose) print(summary(dge_pattern))

    # Extract DEG names and flags
    deg_names <- rownames(dge_profile)[as.logical(dge_pattern)]
    degs_quant <- sum(abs(dge_pattern))

    # Plots log-Fold Change versus log-Concentration.
    if(plot_smear)
    {
        plotSmear(dge_diffs, de.tags=deg_names, pch=pt_chars[2])
        abline(h=c(-1, 1), col="blue")
        title_text <- paste(degs_quant)
        if(!is.null(exp_info) && !is.na(exp_info["Cell"])) title_text <- paste(title_text, exp_info["Cell"])
        title_text <- paste(title_text, "DEGs")
        title_text <- paste(title_text, "in", paste0(group_pair, collapse="/"), "at FDR", fdr)
        if(!is.null(exp_info) && !is.na(exp_info["Group"])) title_text <- paste(title_text, "of", exp_info["Group"], "Group")
        if(!is.null(exp_info) && !is.na(exp_info["Time"])) title_text <- paste(title_text, "at", exp_info["Time"])
        title(main=title_text, cex.main=1.75)
    }

    # Calculate DEGs and save them to files.
    # Generate the combined results of read counts, statistics, and regulation flags for all genes.
    deg_counts_stats <- cbind(dge_counts_stats, dge_pattern)
    col_names <- colnames(deg_counts_stats)
    col_names[length(col_names)] <- "Regulation"
    colnames(deg_counts_stats) <- col_names
    # If deg_only is specified, then only select DEGs for output.
    if(deg_only) deg_counts_stats <- deg_counts_stats[deg_names,]
    # Save DEGs sorted by their range of expression level to data file
    if(!is.null(deg_read_counts_file_name) && deg_read_counts_file_name!="")
    {
        # Generate R-compatible or view-compatible tab-delimited data file.
        if(r_compat) write.table(deg_counts_stats, deg_read_counts_file_name, sep="\t", quote=FALSE)
        else write.table(cbind(Gene=rownames(deg_counts_stats), deg_counts_stats), deg_read_counts_file_name, sep="\t", quote=FALSE, row.names=FALSE)
    }

    # Close plot device
    if(plot_pdf_flag)
    {
        dev.off()
        plot_pdf_flag = FALSE
    }

    # Return the list of DEG read counts.
    return(deg_counts_stats)
}
