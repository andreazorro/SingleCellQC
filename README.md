# Single Cell Quality Control and Integration Pipeline with Seurat
This repository contains an R script for a pipeline designed to automate single-cell quality control and sample integration in Seurat using the output from Cell Ranger. The pipeline utilizes various R packages including dplyr, Seurat, ggplot2, sctransform, glmGamPoi, and tidyverse.

### Usage:
1. **Clone the Repository:**
   ```bash
    git clone https://github.com/andreazorro/SingleCellQC.git
   ```

2. **Place Cell Ranger Output Files:**
   Place your Cell Ranger output folders in the same folder as the R script.

4. **Run the Pipeline:**
   Open your terminal or command prompt, navigate to the folder containing the R script and Cell Ranger output files, and execute the following command:

   ```bash
    Rscript IntegrationDataSeurat.R
   ```

### Output:

Upon successful execution, the pipeline will generate:
  * A folder containing quality control plots.
  * A Seurat object with integrated samples.

## Requirements:
R (version >= 3.5.0)
Required R Packages: 
* dplyr
* Seurat
* ggplot2
* sctransform
* glmGamPoi
* tidyverse

## Notes:
Ensure that all required R packages are installed in your R environment prior to running the pipeline.
Make sure that your Cell Ranger output files are appropriately formatted and contain necessary information for the pipeline to execute successfully.
