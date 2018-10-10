# Convert the read counts and experiment design files by switching
# the position of type component and series component in read counts
# and experiment design file names.
#
# Input file: [main].[series].[type].[ext]
# Output file: [main].[type].[series].[ext]

swapSeriesAndTypeInFileName <- function(name, sep=".")
{
    # Load required library.
    require(tools)

    # Get main and extension name.
    main_name <- file_path_sans_ext(name)
    ext_name <- file_ext(name)

    # Split main name.
    if(sep == ".") main_name_parts <- strsplit(main_name, split="\\.")[[1]]
    else main_name_parts <- strsplit(main_name, split=sep)[[1]]
    n_parts <- length(main_name_parts)

    # If the main name has more than 2 components.
    if(n_parts > 2)
    {
        # Extract name components.
        type_name <- main_name_parts[n_parts]
        series_name <- main_name_parts[n_parts-1]
        # Composite new file name.
        main_name <- paste0(main_name_parts[1:(n_parts-2)], collapse=sep)
        new_name <- paste(main_name, type_name, series_name, ext_name, sep=sep)
    }
    else new_name <- paste(main_name, ext_name, sep=sep)

    # Return the new file name.
    return(new_name)
}
