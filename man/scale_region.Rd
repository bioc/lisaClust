% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hatchingPlot.R
\name{scale_region}
\alias{scale_region}
\alias{scale_region_manual}
\title{Scale constructor for regions}
\usage{
scale_region(aesthetics = "region", ..., guide = "legend")

scale_region_manual(..., values)
}
\arguments{
\item{aesthetics}{The names of the aesthetics that this scale works with}

\item{...}{Arguments passed on to discrete_scale}

\item{guide}{A function used to create a guide or its name. See guides() for more info.}

\item{values}{a set of aesthetic values to map data values to. If this is a
named vector, then the values will be matched based on the names. If unnamed,
values will be matched in order (usually alphabetical) with the limits of the scale.
Any data values that don't match will be given na.value.}
}
\value{
a ggplot guide
}
\description{
Region scale constructor.
}
\examples{

library(spicyR)
## Generate toy data
set.seed(51773)
x <- round(c(
    runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
    runif(200) + 3, runif(200) + 2, runif(200) + 1, runif(200)
), 4) * 100
y <- round(c(
    runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
    runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3
), 4) * 100
cellType <- factor(paste("c", rep(rep(c(1:2), rep(200, 2)), 4), sep = ""))
imageID <- rep(c("s1", "s2"), c(800, 800))
cells <- data.frame(x, y, cellType, imageID)
cells <- SingleCellExperiment::SingleCellExperiment(colData = cells)

## Generate regions
cells <- lisaClust(cells, k = 2)

# Plot the regions with hatchingPlot()
hatchingPlot(cells) +
    scale_region_manual(
        values = c(1, 4), labels = c("Region A", "Region B"),
        name = "Regions"
    )

}
