# Summarize the statistics of experiment design.

sumExperimentDesign <- function(targets)
{
    if(is.null(targets) || (!is.null(targets) && nrow(targets)==0))
    {
        warning("Input experiment design is empty!")
        return(NULL)
    }

    state_cell_plate_rep <- NULL

    # Separate CTRL after drug-treated conditions.
    targets_ctrls <- targets[targets$State=="CTRL",]
    ctrl_cell_plate_names <- paste(targets_ctrls$State, targets_ctrls$Cell, targets_ctrls$Plate, sep=".")
    targets_drugs <- targets[targets$State!="CTRL",]
    drug_cell_plate_names <- paste(targets_drugs$State, targets_drugs$Cell, targets_drugs$Plate, sep=".")

    # 1) Calculate the occurrence of each condition in each cell line.
    # 2) Move CTRL after drug-treated conditions.
    state_cell_plate_rep <- c(table(drug_cell_plate_names), table(ctrl_cell_plate_names))

    # Convert a list of separated condition and cell names to a data frame.
    state_cell_plate_names <- names(state_cell_plate_rep)
    state_cell_plate_rep <- data.frame(t(sapply(strsplit(state_cell_plate_names," *\\. *"),c)), state_cell_plate_rep, row.names=NULL, stringsAsFactors=FALSE)
    colnames(state_cell_plate_rep) <- c("State", "Cell", "Plate", "Replicate")

    # Return the summary.
    return(state_cell_plate_rep)
}
