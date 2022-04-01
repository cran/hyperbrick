#' Spectral Signature of a Hyperspectral Image
#'
#' Visualize statistics calculated through the bands of a hyperspectral image.
#'
#' @param x A numeric matrix or vector containing the values to be plotted at
#' each spectral band (wavelength). Generally, an object obtained with
#' [slideBrick()].
#'
#' @param ... Further graphical parameters. See [par()].
#'
#' @seealso [slideBrick()]
#'
#' @examples
#' p <- system.file('exdata', 'obory.dat', package = 'hyperbrick')
#' im <- buildBrick(p, ref_layer = 35,
#'                 spectral_feature = "radiance",
#'                 hFOV = 36.8, vFOV = 36.8, height = 45)
#' plotRGB(im, r = 63, b = 34, g = 11, scale = 90, axes = TRUE)
#'
#' sw <- slideWindows(im)
#' lapply(sw, lines, col = "white") -> null_obj
#'
#' sb <- slideBrick(im, sw, fun = mean)
#' head(sb)
#'
#' viewSpectra(sb, ylab = "Radiance")
#'
#' @importFrom graphics matplot
#'
#' @aliases viewSpectra
#'
#' @export
viewSpectra <- function(x, ...) {
   if (length(dim(x)) > 1) bands <- rownames(x) else bands <- names(x)
   wavelength <- as.numeric(sub("b", "", bands))
   matplot(wavelength, x, type = "l", ...)
}
