# Biological Analysis Script

This document serves as a guide for the R script used for biological analysis. The script analyzes the effect of various predictor variables on a biological outcome and generates several informative plots.


## Script Execution

1. Ensure you have R and RStudio installed on your system.
2. Open the R markdown in RStudio.
3. Modify the file paths for the data files (replace `../data/male_CS.csv` with the actual path) if needed.
4. Run the markdown file by clicking Knit. It is important to install the knit package as well.

## Output

The script generates several high-resolution PNG images within a `results` folder. Various other plots are also generated. Refer the PDF attached to this repository as well for a much clear explanantion of the analysis.

## Additional Notes

- The script utilizes the `ggplot2` package for creating informative visualizations.
- The script sets specific plot dimensions (width and height) in millimeters (mm) for potential integration into a LaTeX document.
- The script adjusts plot elements like title size, axis labels, and grid lines for better readability.

For further details, refer directly to the R code itself and the R markdown file.
