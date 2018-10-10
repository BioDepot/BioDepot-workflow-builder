# Adapt the field names of read counts datasets to the format used in
# subsequent comparison analysis.

adaptReadCountsDatasets <- function(exprt_design, name_map=NULL)
{
    # name_map: c(Sample="Well",Culture="Experiment",Measure="Plate")

    # Initialization.
    exprt_design_mapped <- exprt_design

    # Change field names according to given name map.
    if(is.vector(name_map) && is.vector(names(name_map)))
    {
        field_names <- names(exprt_design)
        field_names_mapped <- name_map[field_names]
        field_names_mapped_flag <- !is.na(field_names_mapped)
        field_names[field_names_mapped_flag] <- field_names_mapped[field_names_mapped_flag]
        names(exprt_design_mapped) <- field_names
    }

    # Return experiment design table with mapped field names.
    return(exprt_design_mapped)
}
