##### Function to verify and install packages #####
Sys.getenv("R_LIBS_USER")
Sys.getenv("R_LIBS")

if (!require(BiocManager, quietly = TRUE)) {
  install.packages(BiocManager)
}

BiocManager::install(version = "3.20", ask = FALSE)

check_and_install_with_bc <- function(package) {
  if (!require(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

check_and_install_cran <- function(paquete) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    install.packages(paquete, dependencies = TRUE)
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
