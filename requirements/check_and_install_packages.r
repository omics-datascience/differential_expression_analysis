##### Function to verify and install packages #####

# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
#   BiocManager::install(version = "3.20", ask = FALSE)
# }

check_and_install_with_bc <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

check_and_install_cran <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

##### Installing libraries #####

# Packages that are installed using BiocManager
check_and_install_with_bc("limma")
check_and_install_with_bc("edgeR")

# packages that are installed using CRAN repositories
check_and_install_cran("matrixStats")
check_and_install_cran("ggrepel")
check_and_install_cran("dplyr")
check_and_install_cran("ggplot2")
check_and_install_cran("reshape2")
check_and_install_cran("gridExtra")
check_and_install_cran("pROC")

# Load libraries once installed
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(pROC))
