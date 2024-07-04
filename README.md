# Pan Solid Copy Number Variant Validation

This repository is for R scripts used in the ongoing validation of using a Qiagen CLC bioinformatics pipeline (PanSolid 1.6dev) to detect somatic copy number variants in the "Pan Solid" Qiaseq primer panel at the North West Genomic Laboratory Hub (Manchester).

## Sample Identification

The consistent sample identifier is the DNA Database "lab number": an 8 digit number, with the first two digits specifying the year.

Example: 21011525 (the 11,525th sample receieved in 2021).

I have tried to consistently refer to this as "labno", as this is the snakecase version of how it appears on the DNA Database (LABNO).

I have also tried to format labno consistently as a character not a numeric, because this is how it is stored in DNA database and it makes joining tables easier.

## Data

All data is saved internally at the North West GLH, as it contains patient identifiable information. 

**No data should be available in this Github repository.**

The filepath for the data and outputs folders can be set in "set_shared_drive_filepath.R". 
