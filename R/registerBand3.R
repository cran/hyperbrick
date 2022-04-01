#' Single Band-to-Band Registration (Rotation and Translation)
#'
#' Hyperspectral image acquisition normaly causes spatial misalignment between
#' the spectral bands (layers) due to both equipment (such as band-to-band
#' recording delay) and external factors (e.g. sensor vibrations). In this case,
#' a geometric correction is necessary for remote sensing applications such
#' as combining/merging spectral bands. This function uses the HOG (Histogram
#' of Oriented Gradient) descriptor in order to find the optimal rotation
#' angle and translation (xy shift) on a 'slave' band to be spatially align
#' with a 'master' (reference) band.
#'
#' @param slave An object of class \code{RasterLayer} (from package
#' [raster]).
#'
#' @param master An object of class \code{RasterLayer} (from package
#' [raster]).
#'
#' @param ncells An integer giving the number of cells to compute the oriented
#' gradients of the HOG descriptor. Default is 24. See [OpenImageR::HOG()].
#'
#' @param orient An integer giving the number of orientations to compute the
#' oriented gradients of the HOG descriptor. Default is 8. See [OpenImageR::HOG()].
#'
#' @param start_affine A numeric vector containing the starting values for the
#' affine parameters to be optimized, i.e., the shift in x and y and the
#' rotation angle (in degrees). Example: \code{start_affine = c(1, 0, -2)},
#' which indicates a shift of 1 in the x-axis, 0 in the y-axis and a
#' clockwise (negative values) angle of 2 degrees. See more in [affineBrick()].
#'
#' @details This should be used carefully, as rotation affects the spatial
#' dimensions. It is recommended to try [registerBand()] first.
#'
#' The affine parameters are estimated using a general optimization algorithm.
#'
#' @return An object of the same classe as the input \code{slave}, with
#' the fixed extent. An additional attribute called \code{'affine_pars'} is
#' stored, containing the rotation angle (degrees) and the shift in x and y
#' in the same unit as the spatial extent of the image.
#'
#' @seealso [OpenImageR::HOG()], [registerBrick()], [registerBand()]
#'
#' @examples
#' p <- system.file('exdata', 'soybean.tif', package = 'hyperbrick')
#' im <- brick(p)
#' print(im)
#'
#' # see how layer 1 is misregistered
#' plot(im[[3]], col = gray.colors(20), asp = 0)
#' plot(im[[1]], add = TRUE, legend = FALSE,
#'      col = adjustcolor(terrain.colors(20), 0.6))
#'
#' # remove the #s to run
#' # b1_reg <- registerBand3(slave = im[[1]], master = im[[3]],
#' #                         start_affine = c(0, 0, -2.5))
#' # attr(b1_reg, "affine_pars")
#'
#' # plot(im[[3]], col = gray.colors(20), asp = 0)
#' # plot(b1_reg, add = TRUE, legend = FALSE,
#' #      col = adjustcolor(terrain.colors(20), 0.6))
#'
#'
#' @importFrom OpenImageR HOG
#' @importFrom raster extent crop res shift
#' @importFrom stats optim
#'
#' @aliases registerBand3
#'
#' @export
registerBand3 <- function(slave, master,
                          ncells = 24, orient = 8,
                          start_affine)
{
  bx <- res(master)[1] * ncol(master)/10
  by <- res(master)[2] * nrow(master)/10
  ex <- extent(master)[]
  pol <- extent(c(ex[1]+bx, ex[2]-bx, ex[3]+by, ex[4]-by))
  r1c <- raster::as.matrix(crop(master, pol))
  hog1 <- HOG(r1c, cells = ncells,
              orientations = orient)
  fun_affine <- function(par) {
    sx <- par[1]; sy <- par[2]; a <- par[3]
    affpol <- affineCoords(pol, angle = a, c(sx, sy))
    r2c <- raster::as.matrix(crop(slave, extent(affpol)))
    hog2 <- HOG(r2c, cells = ncells,
                orientations = orient)
    log(dist(rbind(hog1, hog2))[1])
  }
  p <- optim(start_affine, fun_affine)
  sx <- p$par[1]; sy <- p$par[2]; a <- p$par[3]
  fixed <- affineBrick(slave, a, c(-sx, -sy))
  attr(fixed, "affine_pars") <- c(sx = sx, sy = sy, angle = a)
  return(fixed)
}
