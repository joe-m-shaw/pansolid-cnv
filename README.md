# Pan Solid Copy Number Variant Validation

This repository is for R scripts used in the ongoing validation of using a Qiagen CLC bioinformatics pipeline to detect somatic copy number variants in the "Pan Solid" Qiaseq primer panel at the North West Genomic Laboratory Hub (GLH) in Manchester.

## Sample Identification

The consistent sample identifier is the DNA Database "lab number": an 8 digit number, with the first two digits specifying the year.

Example: 21011525 (the 11,525th sample receieved in 2021).

I have tried to consistently refer to this as "labno", as this is the snakecase version of how it appears on the DNA Database (LABNO).

I have also tried to format labno consistently as a character not a numeric, because this is how it is stored in DNA database and it makes joining tables easier.

## Data

All data is saved internally at the North West GLH, as it contains patient identifiable information. 

**No data should be available in this Github repository.**

The filepath for the data and outputs folders can be set in "set_shared_drive_filepath.R". 

## Replicating the *ERBB2* validation analysis

To replicate the analysis performed for the validation of detection of *ERBB2* somatic amplifications (summarised in document DOC6260 version 1), follow these steps:

1) Navigate to commit b18d22a in this repo and select "Download ZIP".

2) Save the file contents to your local drive.

3) In the same drive location, copy the data and outputs folders from this shared location: S:/central shared/Genetics/Mol_Shared/Development.Team/Pan Solid CLC Somatic Amplifications Validation/erbb2_validation_r_files

4) Use set_here() to specify the starting point for the here package, as it will default to the scripts folder.

5) In the "erbb2_validation.qmd" file, ammend lines 1361 and 1363 so that "repeatability_fold_change_sd" and "reproducibility_fold_change_sd" become "repeatability_fold_change_sd[[2]]" and "reproducibility_fold_change_sd[[2]]".

6) Render the "erbb2_validation.qmd" script. Check the outputs against the plots and tables in DOC6260 v1. The timestamped
erbb2_validation_dataset_for_export.csv file can be checked against the timestamped data file saved on QPulse using all.equal().

