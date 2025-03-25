# Master's Program in Bioinformatics - Lund University

This repository contains a collection of scripts, data, and reports generated as part of my Master's program in Bioinformatics at Lund University. The content is organized into three subfolders, each corresponding to a specific course focused on coding and computational techniques.

## Subfolders

### 1. BINP16 - Python Programming
This folder contains:
- A variety of Python programming exercises, assignments, and exams.
- Subdirectories for specific running exercises, each with their own:
  - **Data**: Input files such as sequence data, parameter files, and example outputs.
  - **Scripts**: Python scripts used for solving bioinformatics problems, including `malaria.py` and `FastaAligner.py`.
  - **Results**: Output files generated by the scripts, including alignment results and heatmaps.
  - **Visualizations**: Plots and dendrograms generated for comparison studies.

### 2. BIOS13 - Modelling Biological Systems
This folder includes:
- Exam questions and presentation files (PDFs and PowerPoint) for various topics.
- Scripts in R for modeling biological systems, such as population dynamics and gene frequency simulations.
- Exam-related scripts and plots, including phase plane diagrams, time series, and coupled dynamics simulations.

### 3. BIOS14 - Processing and Analysis of Biological Data
This folder encompasses:
- Input data files related to damselfly studies, including CSV files for female and male datasets.
- Scripts in R and related files (RMarkdown, LaTeX, and logs) for statistical analysis and data visualization.
- Results including ANOVA plots, histograms, QQ plots, scatter plots, and contour plots to analyze interactions and thresholds.

### 4. BINP28 - DNA Sequence Informatics I

This repository contains a bash script (`variant_processing.sh`) that performs a comprehensive downstream analysis of a Variant Call Format (VCF) file. The analysis includes variant filtering, population structure inference using ADMIXTURE, and phylogenetic tree construction.  The script utilizes several bioinformatics tools, including `bcftools`, `vcftools`, `plink`, `iqtree2`, `fastme`, and `R`.

This folder includes:
- A scripts folder that contains all the scripts used for the analysis.
- It also contains a README.md file and an environment YML file. 

### 5. BINP29 - Malaria Phylogenetics and Comparative Genomics Project

This folder contains scripts, data, and results related to a bioinformatics pipeline developed for the phylogenetic and comparative genomic analysis of *Plasmodium* species. The core of the project is the `malaria_pipeline.sh` script, which automates genome filtering, gene prediction, ortholog identification, and phylogenetic tree reconstruction using RAxML.

#### **Folder Structure:**
- **Scripts**: Contains essential scripts, including:
  - `malaria_pipeline.sh`: The main pipeline script for processing genomic data.
  - `filtering_scaffolds.py`: Python script for filtering *Haemoproteus tartakovskyi* genome scaffolds based on GC content and length.
  - `busco_protein_fasta.py`: Extracts protein sequences for BUSCO analysis.
  - `datParser.py`: Parses and extracts relevant genomic information from `.dat` files.

#### **Pipeline Overview:**
The `malaria_pipeline.sh` script executes the following steps:
1. **Data Preprocessing:**
   - Calculates genome size and GC content.
   - Extracts gene counts from annotation files.
2. **Genome Filtering:**
   - Removes *Haemoproteus tartakovskyi* scaffolds with GC content >28% and length <3000 nt.
3. **Gene Prediction:**
   - Uses GeneMark-ES for gene prediction.
4. **Functional Annotation:**
   - Performs BLASTx searches against SwissProt to identify homologous genes.
5. **Phylogenetic Analysis:**
   - Aligns sequences and constructs maximum-likelihood phylogenetic trees with RAxML.

A dedicated conda environment (`malaria_env`) manages dependencies, with an additional `busco_env` for BUSCO analysis. The repository ensures reproducibility by specifying software environments and input data locations.

This project provides a structured approach to investigating *Plasmodium* evolution and comparative genomics through automated bioinformatics workflows.

### 6.BINP29 - Ancestry Map - A tool to visualize local acestry

This project focuses on biogeographical analysis by determining an individual's geographical origins based on their DNA at a local ancestry level. Unlike global ancestry, which provides overall genetic composition, local ancestry breaks the genome into smaller chromosomal segments to infer their geographical origins. The Geographical Population Structure (GPS) algorithm, developed by Elhaik et al. (2014), plays a crucial role in this analysis, allowing for a more detailed understanding of genetic heritage.

This pipeline processes Q_files, BIM, and FAM files, runs the GPS algorithm in parallel, merges segment predictions, and generates a dataset for visualization through a web application.

#### Project Workflow

1. **Input Data Preparation**
   - Process Q_files, BIM, and FAM files.
   - Q_files contain admixture proportions for 1069 individuals across 9 gene pools.
   - Output: Prepared dataset for GPS.

2. **Run GPS Algorithm (Parallel Processing)**
   - Process multiple individuals in parallel using the Ray framework.
   - Reduces runtime significantly from hours to ~12 minutes.

3. **Merge Chromosome Predictions**
   - Calculate Haversine distance between adjacent segments.
   - Merge segments based on a user-defined threshold.

4. **Prepare Files for Re-running GPS**
   - Split merged predictions into smaller files (1000 entries each) for efficiency.

5. **Re-run GPS Algorithm on Split Files**
   - Obtain refined geographical coordinates for merged segments.

6. **Final Merge for Plotting Dataset**
   - Adjust segments falling in water bodies.
   - Map each segment to a country and continent.
   - Prepare dataset for visualization.


## Note
The materials in this repository represent only the coding-focused courses I have taken. They include a variety of exercises, assignments, and exams that were submitted and completed during the program. See each folder for a detailed explanation of the projects, exams and exercises.
