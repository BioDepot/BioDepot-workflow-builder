# Record multiple data files generated from DEG comparisons between mulitple pairwise conditions.

recordGeneExpDiffs <- function(deg_read_counts, targets, state_names, read_counts, cell_line, plate, n_top, read_counts_file=NULL, exp_configs_file=NULL, top_stats_file=NULL, deg_stats_file=NULL, sub_char=FALSE, verbose=FALSE, func_dir=NULL)
{
    # Load user-defined functions.
    if(is.null(func_dir)) func_dir <- getwd()
    source(file.path(func_dir, "recordGeneExpDiff.R"), local=TRUE)

    # Initialize top-gene statistics and DEG names.
    top_read_counts <- list()
    deg_states_names_table <- NULL
    top_states_names_table <- NULL

    # Extract DEG statistics from each compared condition pair.
    for(group_pair in names(deg_read_counts))
    {
        # Get the names of two compared conditions.
        group_pair_name <- unlist(strsplit(group_pair, " *\\. *"))
        # Restore the replcaed "+" in the names of list elements to drug names.
        if(sub_char) group_pair_name <- gsub("_", "\\+", group_pair_name)

        # Get the read counts and DEG statistics of two compared conditions.
        deg_read_counts_group_pair <- deg_read_counts[[group_pair]]

        # Prepare the datasets of state 1 and 2.
        # Prepare the dataset of state 1.
        group_pair_state1_idx <- state_names==as.character(group_pair_name)[1]
        targets_state1 <- targets[group_pair_state1_idx,]
        group_names_state1 <- state_names[group_pair_state1_idx]
        read_counts_state1 <- read_counts[,group_pair_state1_idx]
        read_counts_state1 <- read_counts_state1[rowSums(is.na(read_counts_state1))==0,]
        # Prepare the dataset of state 2.
        group_pair_state2_idx <- state_names==as.character(group_pair_name)[2]
        targets_state2 <- targets[group_pair_state2_idx,]
        group_names_state2 <- state_names[group_pair_state2_idx]
        read_counts_state2 <- read_counts[,group_pair_state2_idx]
        read_counts_state2 <- read_counts_state2[rowSums(is.na(read_counts_state2))==0,]

        # Record DEG and TOP statistics from comparisons between two conditions.
        stats_dataset <- recordGeneExpDiff(group_pair=group_pair_name, deg_read_counts=deg_read_counts_group_pair, targets_state1=targets_state1, group_names_state1=group_names_state1, read_counts_state1=read_counts_state1, targets_state2=targets_state2, group_names_state2=group_names_state2, read_counts_state2=read_counts_state2, cell_line=cell_line, plate=plate, n_top=n_top, read_counts_file=read_counts_file, exp_configs_file=exp_configs_file, top_stats_file=top_stats_file, deg_stats_file=deg_stats_file, verbose=verbose)
        top_read_counts[[group_pair]] <- stats_dataset$top_read_counts
        deg_states_names_table <- rbind(deg_states_names_table, stats_dataset$deg_states_names)
        top_states_names_table <- rbind(top_states_names_table, stats_dataset$top_states_names)
    }

    # Return some datasets.
    stats_datasets <- list(top_read_counts=top_read_counts, deg_states_names_table=deg_states_names_table, top_states_names_table=top_states_names_table)
    return(stats_datasets)
}
