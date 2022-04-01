#' Read and Pre-Process Hyperspectral Image
#'
#' Retrieve the raw data from the ENVI file, simultaneously read the
#' header file and build a Brick containing all the layers (spectral bands)
#' of the hyperspectral image. Optionally, do pre-processing steps on row
#' data, such as: radiometric correction using the Dark Object Subtraction
#' (DOS) method, set up the spatial extents using the coordinates of a
#' reference layer and inputs of the camera, convert the row values
#' (digital numbers) to spectral features such as radiance or reflectance.
#'
#' @param path A character giving the path for the ENVI file (.dat).
#' @param path_hdr (Optional, character) The path for the header file (.hdr).
#' @param hFOV (Optional, numeric) The horizontal Field Of View (in degrees) of
#' the camera. Note: this is used to calculate the spatial extent. See Details.
#' @param vFOV (Optional, numeric) The vertical Field Of View (in degrees) of
#' the camera. Note: this is used to calculate the spatial extent. See Details.
#' @param height (Optional, numeric) The flight altitude (in meters). Note:
#' this is be used to calculate the spatial extent. See Details.
#' @param ref_layer (Optional, integer) The reference layer (spectral band)
#' to be used to set up the geographic location of the image.
#' @param spectral_feature A character giving the name of the spectral feature
#' to be loaded. Must be one of the three:
#' \code{"raw"} (digital numbers), \code{"radiance"}, \code{"reflectance"}.
#' See Details.
#' @param reflectance_method A character for selecting the method to calculate
#' the reflectance values. It must be one of the following:
#' \code{"irradiance"} or \code{"white_panel"}. See Details.
#' @param dark_path (Optional) A character giving the path for the ENVI file containing
#' the dark reference raw values, a.k.a. noisy energy, for DOS calibration.
#' See Details.
#' @param dark_quantile A numeric value used to calculate the quantile
#' of the dark reference raw values of each spectral band. Default is
#' 0.25. This will be used for the DOS radiometric correction.
#' @param white_path (Optional) A character giving the path for the ENVI file
#' containing the white panel. See Details.
#'
#' @details The geographical coordinates (xy) registered in the header file
#' are automatically retrieved. If the arguments \code{hFOV}, \code{vFOV} and
#' \code{height} are passed, \code{buildBrick} will automatically compute the
#' UTM zone and the spatial extent, and set up the coordinate reference system.
#'
#' \code{"radiance"} is obtained by multiplying the raw values of each
#' spectral band by the respective "gain" values registered in the header file
#' (if available).
#'
#' If \code{"reflectance"} is passed as value for the argument
#' \code{spectral_feature}, then the value of the next argument,
#' \code{reflectance_method}, will be used to calculate the reflectance
#' values. The \code{"irradiance"} method consists of using the solar
#' irradiance values registered in the header file (if available) at each
#' spectral band as reference for the radiance values reaching the camera
#' sensor. If \code{"white_panel"} is chosen, then the path for the ENVI
#' file containing the image of the white panel must be passed as value for
#' the respective argument. Note that \code{buildBrick} will consider value
#' of the white panel as the maximum.
#'
#' Make sure that the header file (.hdr) has the same name as the
#' ENVI (.dat) file when pointing them at the arguments \code{dark_path}
#' or \code{white_path}.
#'
#' @return A \code{RasterBrick} object (from the package [raster]).
#'
#' @seealso [read_hdr_envi()], [raster::brick()]
#'
#' @examples
#' # Point to an ENVI data
#' path <- system.file('exdata', 'obory.dat', package = 'hyperbrick')
#' # First, let's check the header file
#' # There are 81 bands. First two are noisy.
#' # Irradiance is available.
#' path_hdr <- system.file('exdata', 'obory.hdr', package = 'hyperbrick')
#' readLines(path_hdr)
#'
#' # Example 1 - raw values
#' path <- system.file('exdata', 'obory.dat', package = 'hyperbrick')
#' b <- buildBrick(path)
#' print(b)
#' plot(b, 35)
#'
#' # Example 2 - set up CRS and compute radiance
#' br <- buildBrick(path, hFOV = 36.8, vFOV = 36.8, height = 45,
#'                 ref_layer = 35, spectral_feature = 'radiance')
#' print(br)
#' plot(br, 35)
#'
#' # Example 3 - DOS correction
#' dpath <- system.file('exdata', 'obory_dark.dat', package = 'hyperbrick')
#' brd <- buildBrick(path, hFOV = 36.8, vFOV = 36.8, height = 45,
#'                   ref_layer = 35, spectral_feature = 'radiance',
#'                   dark_path = dpath)
#' print(brd)
#' plot(brd, 35)
#'
#' # Example 4 - compute reflectance
#' bre <- buildBrick(path, hFOV = 36.8, vFOV = 36.8, height = 45,
#'                   ref_layer = 35, spectral_feature = 'reflectance',
#'                   reflectance_method = "irradiance",
#'                   dark_path = dpath)
#' print(bre)
#' plot(bre, 35)
#'
#' # Example 5 - there is a white-reference panel in this image
#' idmax <- which.max(bre[[35]])
#' plot(bre, 35, asp = 0)
#' points(xyFromCell(bre[[35]], idmax), pch = 3)
#'
#' bre2 <- buildBrick(path, hFOV = 36.8, vFOV = 36.8, height = 45,
#'                    ref_layer = 35, spectral_feature = 'reflectance',
#'                    reflectance_method = "white_panel",
#'                    white_path = path, dark_path = dpath)
#' print(bre2)
#' plot(bre2, 35)
#'
#' @importFrom caTools read.ENVI
#' @importFrom stats quantile
#' @importFrom raster brick extent crs
#'
#' @aliases buildBrick
#'
#' @export
buildBrick <- function(path,
                       path_hdr = sub(".dat", ".hdr", path, fixed = TRUE),
                       hFOV = NULL, vFOV = NULL, height = NULL,
                       ref_layer = 1,
                       spectral_feature = c("raw", "radiance", "reflectance"),
                       reflectance_method = c("irradiance", "white_panel"),
                       dark_path = NULL, dark_quantile = 0.25,
                       white_path = NULL)
{
  HDR <- read_hdr_envi(path_hdr, hFOV, vFOV, height)
  dat <- read.ENVI(path, path_hdr)
  A <- array(dat, dim = HDR$dim)
  if (!is.null(dark_path)) {
    dark_raw <- read.ENVI(dark_path,
                          sub(".dat", ".hdr", dark_path,
                              fixed = TRUE))
    dark_cube <- array(dark_raw, dim = HDR$dim)
    dark_vals <- apply(dark_cube, 3, quantile,
                       p = dark_quantile)
    A <- sweep(A, 3, dark_vals)
  }
  feat <- match.arg(spectral_feature)
  if (feat == "radiance") {
    A <- sweep(A, 3, HDR$gain, FUN = "*")
  } else if (feat == "reflectance") {
    refl <- match.arg(reflectance_method)
    if(refl == "irradiance") {
      rad_cube <- sweep(A, 3, HDR$gain, FUN = "*")
      A <- sweep(rad_cube, 3, HDR$irradiance, FUN = "/")
    } else {
      white_raw <- read.ENVI(white_path,
                             sub(".dat", ".hdr", white_path,
                                 fixed = TRUE))
      white_cube <- array(white_raw, dim = HDR$dim)
      if (!is.null(dark_path)) white_cube <- sweep(white_cube, 3, dark_vals)
      white_vals <- apply(white_cube, 3, max)
      A <- sweep(A, 3, white_vals, FUN = "/")
    }
  }
  rb <- brick(A)
  if(!is.null(HDR$extents)) {
    raster::extent(rb) <- raster::extent(HDR$extents[ref_layer,])
    raster::crs(rb) <- HDR$CRS
  }
  names(rb) <- paste0("b", HDR$wavelength)
  return(rb)
}
