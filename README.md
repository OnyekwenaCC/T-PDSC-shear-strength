# T-PDSC-shear-strength
TPDSC model for shear strength &amp; E50 estimation 

**Contents**

TPDSC_shear_strength_E50.R — Main R script

strain_data.xlsx — Experimental strain input data

TPDSC_out_from_excel/ — Output directory (created automatically after running the script)

**How to Run**

Place strain_data.xlsx and the R script in the same directory.

Open R or RStudio and run:
source("TPDSC_shear_strength_E50.R")

Results (peak stresses, E50 values, shear strength parameters) will be exported as CSV files in the output folder.

**License**

This project is released under the MIT License. Please cite appropriately when used in publications.
