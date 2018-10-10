# Generate a gradient from a set of given colors.

generateColorGradient <- function(pal=c("blue","yellow","red"), n=100, test=FALSE)
{
    # n <- 100 # Number of colors to generate

    # Load required library
    require(graphics)

    # Check the validity of specified color palette.
    total_colors <- colors()
    pal_validity <- all(pal %in% total_colors)
    stopifnot(pal_validity)

    # Set color map parameters.
    color_pal <- colorRampPalette(pal)
    color_grad <- color_pal(n)

    # Test the color codes.
    if(test)
    {
        print("Color Values:")
        print(color_grad)
        main <- paste("Color Gradient from", paste0(pal,collapse=", "))
        image(1:n, 1, as.matrix(1:n), col=color_grad, xlab="Gradient", ylab="", main=main, xaxt="n", yaxt="n", bty="n")
    }

    # Return the color set.
    return(color_grad)
}
