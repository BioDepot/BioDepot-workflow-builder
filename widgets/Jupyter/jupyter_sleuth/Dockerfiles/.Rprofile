options(bitmapType="Cairo")
..First <- function() {
  cat("\nSetting default font family to Cairo for base R plots...\n")
  options(device = function(...) {
    Cairo::Cairo(type = "raster", ...)
  })
  par(family = "Cairo")
}
.First <- function() {
  cat("\nSetting default font family to Cairo for ggplot2 plots...\n")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_set(ggplot2::theme_minimal(base_family = "Cairo"))
  }
}
