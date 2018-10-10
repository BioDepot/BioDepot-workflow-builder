# Retrieve a rainbow color set for heat map.

generateRainbowHeatColors <- function(n=100, s=1, v=1, alpha=1, test=FALSE)
{
    # n <- 100 # Number of colors to generate
    # s <- 1 # Contrast between 0 and 1
    # v <- 1 # Brightness between 0 and 1
    # alpha <- 1 # Transparency between 0 and 1
    
    # Load required library    
    require(graphics)
    
    # Generate the rainbow colors for heat map: red for highest and blue for lowest.
    start <- 0
    end <- (max(1, n-1)*2)/(n*3)
    h <- seq.int(start, ifelse(start > end, 1, 0) + end, length.out = n)%%1
    color <- rev(hsv(h, s, v, alpha))
    
    # Test the color codes.
    if(test)
    {
        print("Color Values:")
        print(color)
        image(1:n, 1, as.matrix(1:n), col=color, xlab="Rainbow Heat Colors", ylab="", xaxt="n", yaxt="n", bty="n")
    }
    
    # Return the color set.
    return(color)
}
