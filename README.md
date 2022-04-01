
<img align="right" src="logo_hyperbrick.png" width="16%" height="16%">


*hyperbrick: Accessory Tools for (Pre)Processing Hyperspectral Images*

The package **hyperbrick** is up to read and (pre)process hyperspectral images. These type of sensor data is usually recorded in some raw format. This package contains some easy-to-use functions to promptly build the image with some basic radiometric calibrations and setting up the spatial information. Geometric correction can be done with band-to-band registration (translation and rotation). Further functionalities allows to compute sliding windows statistics over the image.

# Getting started

You can install the development version from GitHub, using:

```r
devtools::install_github("arsilva87/hyperbrick")
```
Afterwards, just load it and it will be ready to use.

```r
library("hyperbrick")
```

Check the package documentation to see examples of its functions.

# Contact

**hyperbrick** is an ongoing project. Contributions are very welcome. If you have a question or have found a bug, please open an [Issue](https://github.com/arsilva87/hyperbrick/issues) or reach out directly by [e-mail](mailto:anderson.silva@ifgoiano.edu.br).
