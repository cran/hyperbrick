# --------------------------------------
# on loading hyperbrick
.onAttach <- function(lib, pkg)
{
  vers <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage(paste("---\nhyperbrick version", vers))
}
