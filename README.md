# Pan Solid Copy Number Variant Validation

This repository is for R scripts used in the ongoing validation of using a Qiagen CLC bioinformatics pipeline to detect somatic copy number variants in the "Pan Solid" Qiaseq primer panel at the North West Genomic Laboratory Hub (GLH) in Manchester.

## Project Structure

### data

All data is saved internally at the North West GLH, as it contains patient identifiable information. 

**No data should be available in this Github repository.**

The filepath for the data and outputs folders can be set in "set_shared_drive_filepath.R". 

### functions

Related functions are grouped together.

### scripts

The scripts folder contains scripts that process, collate and reformat the raw data. Processed data is saved in the relevant folder on the shared drive.

**dev**: the "dev" subfolder stores scripts that are non-essential for performing the final validation analyses. These are scripts that I use for performing one-off tasks or developing new code.

### vignettes

The vignettes folder contains Quarto markdown documents (.qmd) which show the final data analysis.

### shiny

The shiny folder contains the app.R script for auditing the live PanSolid clinical service.

### outputs

Any outputs from the analyses not displayed in the vignettes markdown documents is saved in the "outputs" folder on the shared drive.
