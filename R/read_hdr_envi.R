#' Read Metadata from an ENVI Header File
#'
#' The ENVI image files are accompanied by an ASCII header file (.hdr format)
#' containing the metadata of the image, such as number of samples (rows),
#' lines (columns), number of spectral bands, wavelength, byte order,
#' data type, gain values (if available), irradiance (if available),
#' coordinates etc. This file is necessary to properly read the image data.
#'
#' @param path A character giving the path for the header file (.hdr).
#' @param hFOV (Optional, numeric) The horizontal Field Of View (in degrees) of
#' the sensor. Note: this is used to calculate the spatial extent. See Details.
#' @param vFOV (Optional, numeric) The vertical Field Of View (in degrees) of
#' the sensor. Note: this is used to calculate the spatial extent. See Details.
#' @param height (Optional, numeric) The flight altitude (in meters). Note:
#' this is be used to calculate the spatial extent. See Details.
#'
#' @details The geographical coordinates (xy) registered in the header file
#' are automatically retrieved. Please be aware that \code{read_hdr_envi}
#' uses the string "gps" to search this field in the header file. If necessary,
#' rename this field in the header file.
#'
#' If the arguments \code{hFOV}, \code{vFOV} and
#' \code{height} are passed, \code{read_hdr_envi} will automatically compute
#' the UTM zone and the spatial extent, and set up the coordinate reference
#' system.
#'
#' @return A \code{list} of the following:
#' \describe{
#'   \item{dim}{ An integer vector with the dimensions of the image: number of
#'   columns (x), rows (y) and layers (spectral bands).}
#'   \item{wavelength}{ A numeric vector of the wavelength registered, the same
#'   length as the number of layers (spectral bands).}
#'   \item{gain}{ A numeric vector of gain values of each spectral band.}
#'   \item{irradiance}{ A numeric vector of the solar irradiance values
#'   registered by the corresponding sensor at each spectral band. Please
#'   check the unit (usually W/(m^2 µm sr)) with the manufacturer
#'   specifications.}
#'   \item{coordinates}{ A numeric matrix with the geographic coordinates (xy)
#'   registered by the sensor at each spectral band.}
#'   \item{extents}{ A four-column numeric matrix with the spatial extents
#'   (xmin, xmax, ymin, ymax), in UTM, of each spectral band. This is
#'   \code{NULL} if the optional arguments are not passed.}
#'   \item{CRS}{ The Coordinate Reference System of the spatial extents.
#'   This is \code{NULL} if the optional arguments are not passed.}
#' }
#'
#' @seealso [buildBrick()]
#'
#' @examples
#' # Example 1
#' path_hdr <- system.file('exdata', 'obory.hdr', package = 'hyperbrick')
#' readLines(path_hdr)
#' read_hdr_envi(path_hdr)
#'
#' # Example 2 - set up the CRS to UTM and retrieve extents
#' read_hdr_envi(path_hdr, hFOV = 36.8, vFOV = 36.8, height = 45)
#'
#' @importFrom rgdal project
#' @importFrom utils read.table
#'
#' @aliases read_hdr_envi
#'
#' @export
read_hdr_envi <- function(path,
	hFOV = NULL, vFOV = NULL, height = NULL)
{
   # retrieving spectral attributes
   h <- read.table(path, sep = "=", strip.white = TRUE,
        row.names = NULL, as.is = TRUE, fill = TRUE)
   o <- match(c("samples", "lines", "bands", "wavelength",
      "data gain values", "solar irradiance"), h[,1])
   n_rows <- as.integer(h[o[1], 2])
   n_cols <- as.integer(h[o[2], 2])
   n_bands <- as.integer(h[o[3], 2])
   N <- n_cols * n_rows * n_bands
   wave_length <- as.numeric(gsub("[[:punct:]]", "",
      gsub(".0", "", unlist(strsplit(h[o[4], 2], split = ",")),
         fixed = TRUE)))
   gain <- as.numeric(gsub("}", "",
      gsub("{", "", unlist(strsplit(h[o[5], 2], split = ",")),
         fixed = TRUE), fixed = TRUE))
   irradiance <- as.numeric(gsub("}", "",
      gsub("{", "", unlist(strsplit(h[o[6], 2], split = ",")),
         fixed = TRUE), fixed = TRUE))
   # retrieving coordinates
   gps <- h[grep("gps", h[,1]), 2]
   ss <- unlist(strsplit(gps, "GNRMC,"))[-1]
   locN <- regexpr("A,", ss) + 2
   locE <- regexpr("N,", ss) + 2
   n <- as.numeric(substring(ss, locN, locN + 1)) +
      as.numeric(substring(ss, locN + 2, locN + 3))/60 +
         (as.numeric(substring(ss, locN + 5, locN + 8))/10000)/60
   e <- as.numeric(substring(ss, locE, locE + 2)) +
      as.numeric(substring(ss, locE + 3, locE + 4))/60 +
         (as.numeric(substring(ss, locE + 6, locE + 9))/10000)/60
   coords <- cbind(x = e, y = n)
   # building spatial extents
   if(!is.null(hFOV) & !is.null(vFOV) & !is.null(height)) {
      utm_zone <- floor((coords[1,1] + 180)/6) + 1
      utmproj <- paste0("+proj=utm +zone=", utm_zone,
         " +datum=WGS84 +ellps=intl +units=m +no_defs")
      coords_utm <- project(coords, utmproj)
      xfov <- hFOV*pi/180
      yfov <- vFOV*pi/180
      x_range <- 2*height*tan(xfov/2)
      y_range <- 2*height*tan(yfov/2)
      mat_aux <- kronecker(coords_utm, matrix(c(1, 1), nrow = 1))
      aux_ext <- c(c(0.5, -0.5)*x_range, c(0.5, -0.5)*y_range)
      extents <- sweep(mat_aux, MARGIN = 2, STATS = aux_ext)
      colnames(extents) <- c("xmin", "xmax", "ymin", "ymax")
   } else { extents = NULL; utmproj = NULL}
   # output
   out <- list(dim = c(n_cols, n_rows, n_bands),
      wavelength = wave_length, gain = gain,
      irradiance = irradiance,
      coordinates = as.matrix(coords),
      extents = extents, CRS = utmproj)
   return(out)
}
