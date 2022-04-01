#' Affine Transformation (rotation and shift) on Images
#'
#' Affine transformations are of type \eqn{f(x) = Ax + b}, where \eqn{x} is
#' the spatial coordinates (2D in this case), \eqn{A} is a rotation matrix
#' (it can also include scale/shear parameters, but only rotation is
#' considered here), and \eqn{b} is the translation (xy shift) parameters.
#'
#' @param Brick An object of class \code{RasterBrick}, \code{RasterStack} or
#' \code{RasterLayer} (from package [raster]).
#'
#' @param angle A numeric value of the angle (in degrees, 0-360) to rotate
#' \code{Brick}. A negative value will change the direction of rotation to
#' clockwise.
#'
#' @param xy_shift A numeric vector of length two with the x and y shift (the
#' translation parameters).
#'
#' @return An object of the same class as the input \code{Brick}.
#'
#' @note Affine transformation affects the image dimension.
#'
#' @seealso [affineCoords()], [registerBrick()]
#'
#' @examples
#' p <- system.file('exdata', 'soybean.tif', package = 'hyperbrick')
#' im <- brick(p)
#' print(im)
#'
#' # view band-3
#' plot(im[[3]], col = gray.colors(20), asp = 0)
#'
#' # rotate band-3 at 3.5 degrees counter-clockwise
#' b3_rot <- affineBrick(im[[3]], angle = 3.5)
#' plot(b3_rot, add = TRUE, legend = FALSE,
#'     col = adjustcolor(terrain.colors(20), 0.5))
#'
#' @importFrom raster brick extent crs rasterize shift values
#'
#' @aliases affineBrick
#'
#' @export
affineBrick <- function(Brick, angle = 0, xy_shift = c(0, 0))
{
   stopifnot(angle >= -360 & angle <= 360)
   if (angle != 0 & angle != 360 & angle != -360) {
      rxy <- affineCoords(Brick, angle, xy_shift)
      rb <- brick()
      raster::extent(rb) <- raster::extent(rxy)
      raster::crs(rb) <- raster::crs(Brick)
      newb <- rasterize(x = rxy, y = rb,
                        field = values(Brick))
   } else {
      newb <- shift(Brick, xy_shift[1], xy_shift[2])
   }
   return(newb)
}
