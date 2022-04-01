#' Affine Transformation (rotation and shift) on Spatial Coordinates
#'
#' Affine transformations are of type \eqn{f(x) = Ax + b}, where \eqn{x} is
#' the spatial coordinates (2D in this case), \eqn{A} is a rotation matrix
#' (it can also include scale/shear parameters, but only rotation is
#' considered here), and \eqn{b} is the translation (xy shift) parameters.
#'
#' @param s A spatial object where [raster::coordinates()] can be extracted from.
#' Example classes: \code{Extent}, \code{RasterLayer}, \code{RasterBrick},
#' \code{RasterStack} (from package [raster]), \code{SpatialPolygon}
#' (from package sp).
#'
#' @param angle A numeric value of the angle (in degrees, 0-360) to rotate
#' \code{s}. A negative value will change the direction of rotation to
#' clockwise.
#'
#' @param xy_shift A numeric vector of length two with the x and y shift (the
#' translation parameters).
#'
#' @return A two-column matrix with the transformed coordinates (xy).
#'
#' @seealso [affineBrick()]
#'
#' @examples
#' p <- system.file('exdata', 'soybean.tif', package = 'hyperbrick')
#' im <- brick(p)
#' print(im)
#'
#' # view band-3
#' plot(im[[3]], col = gray.colors(20), asp = 0)
#'
#' # draw a spatial polygon on image
#' pol <- Polygon(extent(c(40, 85, 50, 150)))
#' lines(pol)
#'
#' # rotate and shift the spatial polygon
#' new_pol <- affineCoords(pol, angle = -3, xy_shift = c(-11, 0))
#' plot(im[[3]], col = gray.colors(20), asp = 0)
#' lines(new_pol)
#'
#' # do some analysis within it, like:
#' new_pols <- SpatialPolygons(list(Polygons(list(Polygon(new_pol)), "id0")))
#' plot(mask(im[[3]], new_pols))
#' mean(extract(im[[3]], new_pols)[[1]])
#'
#' @importFrom raster coordinates
#'
#' @aliases affineCoords
#'
#' @export
affineCoords <- function(s, angle = 0, xy_shift = c(0, 0))
{
   stopifnot(angle >= -360 & angle <= 360)
   xy <- coordinates(s)
   a <- angle*pi/180
   cen <- colMeans(xy)
   aff_mat <- matrix(c(cos(a), -sin(a), cen[1] + xy_shift[1],
      sin(a), cos(a), cen[2] + xy_shift[2]), nrow = 3)
   xy_cen <- sweep(xy, 2, cen)
   new_xy <- cbind(xy_cen, 1) %*% aff_mat
   return(new_xy)
}
