# Calculate differentially expressed molecules from their read counts data,
# such as genes and proteins.

calcExpressionDifference <- function(read_counts, exprt_design, drug_groups, bio_spec, time_length, results_dir, fdrs=c(0.05), state_name_ctrl="CTRL", min_samples=3, rep_req="all", dispersion="auto", n_top=100, deg_only=FALSE, norm_base=1e+6, remove_small=TRUE, plot_bcv=TRUE, plot_smear=TRUE, pt_chars=c(16,19), sub_char=TRUE, read_counts_title="ReadCounts", deg_title="DEG", top_title="TOP", exp_configs_title="ExpConfigs", calc_title="Calc", plot_title="Plots", verbose1=FALSE, verbose2=FALSE, func_dir=NULL)
{
    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "compareGeneExpDiffs.R"), local=TRUE)
    source(file.path(func_dir, "recordGeneExpDiffs.R"), local=TRUE)
    source(file.path(func_dir, "sumDEGCounts.R"), local=TRUE)

    # FDRs for calculating DEGs: 0.05 and 0.1.
    deg_counts_sums <- NULL
    for(fdr in fdrs)
    {
        if(verbose1) print(paste("FDR:", fdr))

        # Data directory for the results at current FDR.
        fdr_dir_name <- paste("FDR", fdr, sep="-")
        fdr_dir <- file.path(results_dir, fdr_dir_name)
        write_file <- !all(c(is.null(read_counts_title), is.null(deg_title), is.null(top_title), is.null(exp_configs_title), is.null(calc_title), is.null(plot_title)))
        if(write_file && !dir.exists(fdr_dir)) dir.create(path=fdr_dir, mode="0755")

        # File names of DEG and TOP statistics.
        if(!is.null(deg_title)) deg_states_names_file <- file.path(fdr_dir, paste(paste(deg_title, fdr, sep="-"), "tsv", sep="."))
        else deg_states_names_file <- NULL
        if(!is.null(top_title)) top_states_names_file <- file.path(fdr_dir, paste(paste(top_title, n_top, sep="-"), "tsv", sep="."))
        else top_states_names_file <- NULL

        # Initialize a list of DEG names for all drug groups.
        deg_read_counts <- list()
        deg_states_names_table <- NULL
        top_read_counts <- list()
        top_states_names_table <- NULL

        # Identify the DEGs by control-drug comparison for all drug groups.
        for(drug_group_name in names(drug_groups))
        {
            # Skip CTRL group if needed.
            if(drug_group_name == state_name_ctrl) next

            # Prepare experimental information.
            exp_info <- c(group=drug_group_name)

            if(verbose1) print(paste(drug_group_name, "Group"))
            drug_group_idx <- rep(FALSE, length(exprt_design$State))
            drug_group <- drug_groups[[drug_group_name]]
            for(drug_item in drug_group) drug_group_idx <- (drug_group_idx | exprt_design$State==drug_item)

            # Initialize a list of DEG names for current drug group.
            deg_read_counts_group <- list()
            deg_states_names_group_table <- NULL
            top_read_counts_group <- list()
            top_states_names_group_table <- NULL

            #bio_spec <- unique(exprt_design$Species)[1]
            #time_length <- paste("Hour", unique(targets_merged$Time)[1], sep=".")
            cell_lines <- levels(factor(exprt_design$Cell))
            for(cell_line in cell_lines)
            {
                cell_type <- paste(bio_spec, cell_line, sep=".")
                if(verbose1) print(paste("In", cell_type))

                # Prepare experimental information.
                exp_info <- c(exp_info, cell=cell_type)
                cell_line_idx <- (exprt_design$Cell==cell_line)
                ctrl_idx <- (exprt_design$State==state_name_ctrl)

                # Initialize a list of DEG names for current cell line.
                deg_read_counts_cell <- list()
                deg_states_names_cell_table <- NULL
                top_read_counts_cell <- list()
                top_states_names_cell_table <- NULL

                plates <- levels(factor(exprt_design$Plate))
                for(plate in plates)
                {
                    plate_code <- paste0("P", plate)
                    plate_name <- paste("Plate", plate, sep=".")
                    plate_idx <- (exprt_design$Plate==plate)

                    ################### Prepare comparison of CTRLs and drugs in current drug group ###################

                    # 1. Prepare targets.
                    # Prepare the dataset of CTRL.
                    targets_ctrl_cell_plate <- exprt_design[ctrl_idx&cell_line_idx&plate_idx,]
                    # Prepare the dataset of all drugs in current drug group.
                    targets_drug_group_cell_plate <- exprt_design[drug_group_idx&cell_line_idx&plate_idx,]
                    if(nrow(targets_ctrl_cell_plate)<min_samples || nrow(targets_drug_group_cell_plate)<min_samples)
                    {
                        # At least min_samples CTRL sample is needed for comparison in a cell line.
                        if(nrow(targets_ctrl_cell_plate)<min_samples) warning(paste("At least min_samples", state_name_ctrl, "sample is required for cell line", paste0(cell_type, "!")))
                        # At least min_samples number of drug samples are needed for comparison in a cell line.
                        if(nrow(targets_drug_group_cell_plate)<min_samples) warning(paste("At least", min_samples, "drug samples are required for cell line", paste0(cell_type, "!")))
                        next
                    }
                    # Generate the targets containing CTRL and current drug group.
                    targets_compare_group_cell_plate <- rbind(targets_ctrl_cell_plate, targets_drug_group_cell_plate)

                    # 2. Prepare state names and read counts.
                    # Generate the state names for comparing CTRL and current drug group.
                    state_names_compare <- targets_compare_group_cell_plate$State
                    # Get the state names of current drug group.
                    state_names_drug <- state_names_compare[state_names_compare!=state_name_ctrl]
                    # Get the read counts matrix of CTRL and current drug group.
                    read_counts_group <- read_counts[,targets_compare_group_cell_plate$ID]
                    if(is.vector(read_counts_group)) read_counts_group <- as.matrix(read_counts_group)
                    colnames(read_counts_group) <- targets_compare_group_cell_plate$ID
                    # Generate the pairs of state names to compare.
                    state_names_drug_unique <- unique(state_names_drug)
                    compared_pair <- data.frame(Group1=rep(state_name_ctrl,length(state_names_drug_unique)), Group2=state_names_drug_unique, stringsAsFactors=FALSE)

                    ################### Set up output file names ###################

                    # File names for read counts and DEG plots.
                    cell_time_plate_title <- paste(cell_type, time_length, plate_name, sep="-")
                    # File names of calculation results containing read counts, normalized
                    # read counts, read counts statistics, and DEG statistics.
                    if(!is.null(calc_title)) calc_file <- file.path(fdr_dir, paste(paste(cell_time_plate_title, calc_title, sep="-"), "tsv", sep="."))
                    else calc_file <- NULL
                    # File names of DEG plots.
                    if(!is.null(plot_title)) plots_file <- file.path(fdr_dir, paste(paste(cell_time_plate_title, plot_title, sep="-"), "pdf", sep="."))
                    else plots_file <- NULL

                    ################### Compare the difference between two conditions ###################

                    # Calculate the DEGs in current experiment.
                    deg_read_counts_plate <- compareGeneExpDiffs(read_counts=read_counts_group, group_names=state_names_compare, group_pairs=compared_pair, dispersion=dispersion, min_rep=min_samples, rep_req=rep_req, fdr=fdr, deg_only=deg_only, norm_base=norm_base, remove_small=remove_small, exp_info=exp_info, deg_read_counts_file=calc_file, plot_bcv=plot_bcv, plot_smear=plot_smear, deg_plots_file=plots_file, pt_chars=pt_chars, sub_char=sub_char, verbose=(if(verbose1) verbose2 else verbose1), func_dir=func_dir)
                    # Save the DEG statistics of current plate of current cell line of current drug group.
                    deg_read_counts_cell[[plate_code]] <- deg_read_counts_plate

                    ################### Generate DEG and TOP statistics ###################

                    # Prepare file names.
                    if(!is.null(read_counts_title)) read_counts_file <- file.path(fdr_dir, paste(paste(cell_time_plate_title, read_counts_title, sep="-"), "tsv", sep="."))
                    else read_counts_file <- NULL
                    if(!is.null(exp_configs_title)) exp_configs_file <- file.path(fdr_dir, paste(paste(cell_time_plate_title, exp_configs_title, sep="-"), "tsv", sep="."))
                    else exp_configs_file <- NULL
                    if(!is.null(top_title)) top_stats_file <- file.path(fdr_dir, paste(paste(cell_time_plate_title, paste(top_title, n_top, sep="."), sep="-"), "tsv", sep="."))
                    else top_stats_file <- NULL
                    if(!is.null(deg_title)) deg_stats_file <- file.path(fdr_dir, paste(paste(cell_time_plate_title, paste(deg_title, fdr, sep="."), sep="-"), "tsv", sep="."))
                    else deg_stats_file <- NULL

                    # Record DEG and TOP statistics from comparisons between mulitple conditions.
                    stats_datasets <- recordGeneExpDiffs(deg_read_counts=deg_read_counts_plate, targets=targets_compare_group_cell_plate, state_names=state_names_compare, read_counts=read_counts_group, cell_line=cell_line, plate=plate, n_top=n_top, read_counts_file=read_counts_file, exp_configs_file=exp_configs_file, top_stats_file=top_stats_file, deg_stats_file=deg_stats_file, sub_char=sub_char, verbose=verbose1, func_dir=func_dir)
                    # Save the read counts statistics of top genes of current plate of current cell line of current group.
                    top_read_counts_cell[[plate_code]] <- stats_datasets$top_read_counts
                    # Save the DEG names of current plate of current cell line of current group.
                    deg_states_names_plate_table <- stats_datasets$deg_states_names_table
                    deg_states_names_cell_table <- rbind(deg_states_names_cell_table, deg_states_names_plate_table)
                    # Save the TOP names of current plate of current cell line of current group.
                    top_states_names_cell_table <- rbind(top_states_names_cell_table, stats_datasets$top_states_names_table)

                    ################### Print DEG and TOP summary ###################

                    # Print the number of unique DEGs found for current plate of current cell line of current group.
                    if(verbose1) print(table(paste(deg_states_names_plate_table$Condition1, deg_states_names_plate_table$Condition2, sep=".vs.")))
                    if(verbose1) print(paste("Totally", length(unique(deg_states_names_plate_table$Gene)), "unique DEGs are found in comparing", paste0(state_names_drug_unique, collapse=" and "), "of", cell_type))
                }

                # Save the DEG statistics of current cell line of current drug group.
                deg_read_counts_group[[cell_line]] <- deg_read_counts_cell
                deg_states_names_group_table <- rbind(deg_states_names_group_table, deg_states_names_cell_table)
                # Save the TOP names and statistics of current cell line of current drug group.
                top_read_counts_group[[cell_line]] <- top_read_counts_cell
                top_states_names_group_table <- rbind(top_states_names_group_table, top_states_names_cell_table)
            }

            # Save the DEG names and statistics of current drug group.
            deg_read_counts[[drug_group_name]] <- deg_read_counts_group
            deg_states_names_table <- rbind(deg_states_names_table, deg_states_names_group_table)
            # Save the TOP names and statistics of current drug group.
            top_read_counts[[drug_group_name]] <- top_read_counts_group
            top_states_names_table <- rbind(top_states_names_table, top_states_names_group_table)
        }

        # Save the DEG names of current drug group.
        if(!is.null(deg_states_names_file)) write.table(deg_states_names_table, deg_states_names_file, sep="\t", quote=FALSE, row.names=FALSE)
        # Save the TOP names of current drug group.
        if(!is.null(top_states_names_file)) write.table(top_states_names_table, top_states_names_file, sep="\t", quote=FALSE, row.names=FALSE)

        #################### Summarize DEG statistics ####################

        if(!is.null(deg_title)) deg_counts_file <- file.path(fdr_dir, paste(paste(deg_title, fdr, "Summary", sep="-"), "tsv", sep="."))
        else deg_counts_file <- NULL
        deg_counts_sum <- sumDEGCounts(deg_stats=deg_states_names_table, deg_counts_file=deg_counts_file, verbose=verbose1, func_dir=func_dir)
        # Summarize DEG statistics at all FDRs.
        deg_counts_sums <- rbind(deg_counts_sums, data.frame(FDR=rep(fdr,nrow(deg_counts_sum)), deg_counts_sum, stringsAsFactors=FALSE))
    }

    # Save DEG statistics at all FDRs.
    if(!is.null(deg_title))
    {
        deg_counts_sums_file <- file.path(results_dir, paste(paste(deg_title, "Summary", sep="-"), "tsv", sep="."))
        write.table(deg_counts_sums, deg_counts_sums_file, sep="\t", quote=FALSE, row.names=FALSE)
    }

    # Return the list of raw read counts and differential comparison statistics of all genes.
    return(deg_read_counts)
}
