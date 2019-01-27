[![DOI](https://zenodo.org/badge/100089287.svg)](https://zenodo.org/badge/latestdoi/100089287)

# SCBursts - Single Channel Bursts Analysis

This `R` Package was designed for use by the [daCosta lab](http://www.dacosta.net/). The package contains the following features:

- Transform data between formats (`.evt`, `.dwt`, etc)
- Correct for effects of the Gaussian filter on the time-series data.
- Identify bursts in time series, and isolate those bursts.
- Provide toolkit for sorting, filtering, and studying individual bursts.
- **In the future** Burst detection.

If any features seem wrongfully absent, or if any methods can be improved upon, feel free to either create an [issue](https://github.com/dacostalab/scbursts/issues) or (even better) submit a pull request.

# Installation

## From a tarball (easy)

Get the `tar.gz` file from [the daCosta lab](http://www.dacosta.net/publications.html). You should then be able to install from this file within R-Commander, R-Studio, or from the command-line with

~~~
R CMD INSTALL scbursts_*.*.tar.gz
~~~

## From CRAN (easier)

Once the package gets uploaded to CRAN, in `R` you can run

```
install.packages("scbursts")
```

and then, as usual, just run

```
library(scbursts)
```

to get started.

## From `R` (with pdf documentation)

If you have `LaTeX` and `pandoc` installed, you can install the package **and** generate some documentation for it (this is optional. If you don't have those see the next section). 

Open an R console, and run the following lines

```{.R}
install.packages("devtools")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("roxygen2")
library(devtools)
install_github("dacostalab/scbursts", build_vignettes = TRUE)

# to update install run
# install_github("dacostalab/scbursts", build_vignettes = TRUE, force=TRUE)
```

You should then be able to call

```{.R}
library(scbursts)
```

and to see documentation run

```{.R}
vignette("scbursts")
```

## From `R`

If you have `LaTeX` and `pandoc` installed, you can install the package **and** generate some documentation for it (this is optional. If you don't have those see the next section). 

Open an R console, and run the following lines

```{.R}
install.packages("devtools")
library(devtools)
install_github("dacostalab/scbursts", build_vignettes = TRUE)

# to update install run
# install_github("dacostalab/scbursts", force=TRUE)
```

You should then be able to call

```{.R}
library(scbursts)
```

## From Source 

**You only want to do this if you plan to modify the package. Otherwise use another option**

1. Start by installing the dependencies. On Ubuntu

```
# apt-get install texlive-full pandoc pandoc-citeproc make r-base pkg-config libcurl4-openssl-dev libxml2-dev
```

With those installed, you can install this from source with `make`. You will need to make sure that you have this installed. Once you have it, the steps are:

2. Get a copy of this repository, either by downloading a zip or by `git clone`-ing.

3. Open a terminal in the directory

4. You will also need some `R` packages. If you have write access to your R library path, then you can run `make deps`. Otherwise in an `R` terminal run

```{.R}
install.packages("devtools")
install.packages("knitr")
install.packages("rmarkdown")
install.packages("roxygen2")
install.packages("gdata")
```

5. Now, back in the shell, you can run `make` and `make install`.

And then the package should be installed. Note that you won't be able to create any manuals unless you have LaTeX and pandoc installed.

# Manual

You can view the soft-documentation for this package by calling

```{.R}
vignette("scbursts")
```

from an R console.
