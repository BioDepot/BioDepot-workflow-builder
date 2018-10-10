# Generate a table of experiment design fields from the header line of
# a read-counts data file, according to a set of standard fields for
# experiment design table, including:
#    Species:   biological species of a sample column.
#    Subject:   species subject of a sample column.
#    State:     biological or perturbation condition of a sample column.
#    Culture:   culture code of a sample column.
#    Replicate: replicate code of a sample column.
#    Measure:   measure code of a sample column.
#    Sample:    sample code of a sample column.
#    Time:      time point of a sample column measurement

generateExperimentDesignTable <- function(read_counts_column_names, split_char=NULL, field_names=NULL, field_types=c(Species="character", Subject="character", State="character", Culture="character", Replicate="character", Measure="character", Sample="character", Time="character"), species_default="Human", subject_default="A", state_default="CTRL", culture_default=1, replicate_default=1, measure_default=1, sample_default=1:length(read_counts_column_names), time_default=48)
{
    # Verify the given field names is contained in standard experiment design fields.
    exprt_design_field_names <- c("Species", "Subject", "State", "Culture", "Replicate", "Measure", "Sample", "Time")
    if(!is.null(split_char))
    {
        stopifnot(!is.null(field_names))
        stopifnot(all(field_names %in% exprt_design_field_names))
        stopifnot(all(names(field_types) %in% exprt_design_field_names))
    }

    # read_counts_column_names: the name of multiple read-counts columns in a read-counts file.
    exprt_design_table <- NULL
    for(col_idx in 1:length(read_counts_column_names))
    {
        # Initialize an array of epxeriment design fields.
        exprt_design_fields <- c(Species=species_default, Subject=subject_default, State=state_default, Culture=culture_default, Replicate=replicate_default, Measure=measure_default, Sample=sample_default[col_idx], Time=time_default)

        # Split the name of a sample coumn into multiple parts to be assigned
        # to standard data fields.
        if(!is.null(split_char))
        {
            # Split the name of a read-counts column into multiple parts.
            field_parts <- strsplit(read_counts_column_names[col_idx], split=split_char)[[1]]
            n_fields <- length(field_names)
            # Verify the number of splitted field parts is equal to that of given field names.
            stopifnot(length(field_parts) == n_fields)
            # Assign splitted field parts to experiment design fields.
            for(field_idx in 1:n_fields) exprt_design_fields[field_names[field_idx]] <- field_parts[field_idx]
        }
        else
        {
            # If split_char.is NULL, each column name won't be splitted and therefore
            # will be used as Sample field solely.
            exprt_design_fields["Sample"] <- read_counts_column_names[col_idx]
        }
        # Add experiment design fields of current read-counts column to experiment design table.
        exprt_design_table <- rbind(exprt_design_table, exprt_design_fields)
    }
    # Convert experiment design table to data frame type.
    exprt_design_table <- data.frame(exprt_design_table, stringsAsFactors=FALSE, row.names=1:nrow(exprt_design_table))

    # Convert each field of experiment design table to specified data type.
    if(!is.null(field_types))
    {
        for(field_name in names(exprt_design_table))
        {
            field_type <- field_types[field_name]
            if(field_type == "character") exprt_design_table[,field_name] <- as.character(exprt_design_table[,field_name])
            else if(field_type == "numeric") exprt_design_table[,field_name] <- as.numeric(exprt_design_table[,field_name])
            else if(field_type == "logical") exprt_design_table[,field_name] <- as.logical(exprt_design_table[,field_name])
            else warning(paste(field_type, "isn't a supported field type!"))
        }
    }
    else warning("Skip field type conversion as field_type isn't specified!")

    # Return generated experiment design table.
    return(exprt_design_table)
}
