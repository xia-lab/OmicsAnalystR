
# Description

**_OmicsAnalystR_** is the underlying R package synchronized with OmicsAnalystR web server. It is designed for data-driven multi-omics integration and systems-level interpretation. The R package is composed of R functions necessary for the web-server to perform data processing, batch correction, correlation analysis, clustering analysis and dimension reduction analysis.

Following installation and loading of _OmicsAnalystR_, users will be able to reproduce web server results from their local computers using the R command history downloaded from OmicsAnalyst. Running the R functions will allow more flexibility and reproducibility.

Note - OmicsAnalystR is still under development - we cannot guarantee full functionality
# Installation

**Step 1. Install package dependencies**

To use OmicsAnalystR, make sure your R version is >4.0.3 and install all package dependencies. Ensure that you are able to download packages from Bioconductor. To install package dependencies, use the pacman R package. Note that some of these packages may require additional library dependencies that need to be installed prior to their own successful installation.

```
install.packages("pacman")

library(pacman)

pacman::p_load(igraph, RColorBrewer, qs, rjson, RSQLite)
```

**Step 2. Install the package**

OmicsAnalystR is freely available from GitHub. The package documentation, including the vignettes for each module and user manual is available within the downloaded R package file. If all package dependencies were installed, you will be able to install the OmicsAnalystR. 

Install the package directly from github using the _devtools_ package. Open R and enter:

```

# Step 1: Install devtools
install.packages(devtools)
library(devtools)

# Step 2: Install OmicsAnalystR WITHOUT documentation
devtools::install_github("xia-lab/OmicsAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))

# Step 2: Install OmicsAnalystR WITH documentation
devtools::install_github("xia-lab/OmicsAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```

# Tips for using the OmicsAnalystR package

1. The first function that you will use in every module is the `Init.Data` function, which initiates the _dataSet_ object that stores user's data for further processing and analysis.
2. The OmicsAnalystR package will output data files/tables/analysis/networks outputs in your current working directory.
3. Every function must be executed in sequence as it is shown on the R Command history, please do not skip any commands as this can result in errors downstream.
4. Each main function in OmicsAnalystR is documented. Use the _?Function_ format to open its documentation. For instance, use `?OmicsAnalystR::QueryNet` to find out more about this function.

# Example

## Processing and integrating datasets

```
library(OmicsAnalystR)

# Step 1. Initiate the R objects
Init.Data();

# Step 2. Read datasets
ReadOmicsData("preg_prot.csv", "prot");
ReadOmicsData("preg_met.csv", "met_t");

# Step 3. Process datasets individually
SanityCheckData("preg_prot.csv");
CheckDataType("preg_prot.csv", "true");
AnnotateGeneData("preg_prot.csv", "hsa", "symbol");
RemoveMissingPercent("preg_prot.csv", 0.5)
ImputeMissingVar("preg_prot.csv", "min")
FilteringData("preg_prot.csv","pct","2", "15");
NormalizingData("preg_prot.csv", "log", "NA", "AutoNorm");

SanityCheckData("preg_met.csv");
CheckDataType("preg_met.csv", "true");
AnnotateMetaboliteData("preg_met.csv", "name");
RemoveMissingPercent("preg_met.csv", 0.5)
ImputeMissingVar("preg_met.csv", "min")
FilteringData("preg_met.csv","pct","2", "15");
NormalizingData("preg_met.csv", "log", "NA", "AutoNorm");

# Step 4. Visual inspection of processed data
PlotMultiTsne("qc_multi_tsne_0_","72", "png", "");
PlotMultiPCA("qc_multi_pca_0_","72", "png", "");
PlotMultiDensity("qc_multi_density_0_","72", "png", "");

# Step 5. Network correlation analysis
DoFeatSelectionForCorr("default",20,3);
DoOmicsCorrelation("univariate", "pearson");
DoCorrelationFilter("both", "true", "both",0.9,0.5,2000.0,"genus","agora","false");
```


