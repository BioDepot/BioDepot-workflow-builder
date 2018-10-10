# Generate drug-group table from an experiment design file and a drug property file.

generateDrugGroupTable <- function(exprt_design_file, drug_props_file, default_group="UnkownGroup", default_type="UnkownType", drug_group_file=NULL)
{
    # Read Experiment design.
    stopifnot(file.exists(exprt_design_file))
    exprt_design_table <- read.delim(file=exprt_design_file, quote="", comment.char="#", stringsAsFactors=FALSE)

    # Read drug properties.
    stopifnot(file.exists(drug_props_file))
    drug_props_table <- read.delim(file=drug_props_file, quote="", comment.char="#", stringsAsFactors=FALSE)
    stopifnot(all(!duplicated(drug_props_table$Drug)))
    rownames(drug_props_table) <- drug_props_table$Drug
    drug_props_table$Drug <- NULL

    # Retrieve the properties of the drugs in experiment design.
    unique_states <- unique(exprt_design$State)
    unique_states <- unique_states[unique_states!=state_name_ctrl]
    drug_group_type_table <- drug_props_table[unique_states, c("Group","Type")]
    drug_group_type_table$Group[is.na(drug_group_type_table$Group)] <- default_group
    drug_group_type_table$Type[is.na(drug_group_type_table$Type)] <- default_type

    # Create a drug-group table.
    drug_group_table <- data.frame(Drug=unique_states, Group=drug_group_type_table$Group, Type=drug_group_type_table$Type, stringsAsFactors=FALSE)
    if(!is.null(drug_group_file)) write.table(drug_group_table, file=drug_group_file, quote=FALSE, sep="\t", na="", row.names=FALSE)

    # Return the drug-group table.
    return(drug_group_table)
}
