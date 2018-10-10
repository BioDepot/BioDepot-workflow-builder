# Convert a experiment design table to a set of standard experiment design
# fields, according to:
#    Species:   biological species of a sample column.
#    Subject:   species subject of a sample column.
#    State:     biological or perturbation condition of a sample column.
#    Culture:   culture code of a sample column.
#    Replicate: replicate code of a sample column.
#    Measure:   measure code of a sample column.
#    Sample:    sample code of a sample column.
#    Time:      time point of a sample column measurement

convertExperimentDesignTable <- function(exprt_design, field_names_positions, field_types=c(Species="character", Subject="character", State="character", Culture="character", Replicate="character", Measure="character", Sample="character", Time="character"), species_default="Human", subject_default="A", state_default="CTRL", culture_default=1, replicate_default=1, measure_default=1, sample_default=1:nrow(exprt_design), time_default=48)
{
    # field_names_positions includes column indices of specified experiment
    # desgin fields in given experiment design table, in a form of:
    # c(Species=1, Subject=2, State=3, Culture=4, Replicate=5, Measure=6, Sample=7, Time=8)

    # Verify the given field names is contained in standard experiment design fields.
    exprt_design_field_names <- c("Species", "Subject", "State", "Culture", "Replicate", "Measure", "Sample", "Time")
    field_names <- names(field_names_positions)
    stopifnot(all(field_names %in% exprt_design_field_names))
    stopifnot(all(names(field_types) %in% exprt_design_field_names))

    # Generate a experiment design table with default values.
    n_samples <- nrow(exprt_design)
    exprt_design_table <- data.frame(Species=rep(species_default,n_samples), Subject=rep(subject_default,n_samples), State=rep(state_default,n_samples), Culture=rep(culture_default,n_samples), Replicate=rep(replicate_default,n_samples), Measure=rep(measure_default,n_samples), Sample=sample_default, Time=rep(time_default,n_samples), stringsAsFactors=FALSE)

    # Assign specified fields of given experiment design to the table.
    for(field_idx in 1:length(field_names_positions))  exprt_design_table[,field_names[field_idx]] <- exprt_design[,field_names_positions[field_idx]]

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

    # Return converted experiment design table.
    return(exprt_design_table)
}
