# Downstream Analysis of VCF File - Population Structure Analysis

### Author
-[@ChandrashekarCR](https://github.com/ChandrashekarCR)

This repository contains a bash script (`variant_processing.sh`) that performs a comprehensive downstream analysis of a Variant Call Format (VCF) file. The analysis includes variant filtering, population structure inference using ADMIXTURE, and phylogenetic tree construction.  The script utilizes several bioinformatics tools, including `bcftools`, `vcftools`, `plink`, `iqtree2`, `fastme`, and `R`.

## Computational Environment

A dedicated conda environment, `vcf_env`, is used to manage the software dependencies.  This ensures reproducibility and avoids conflicts with other software installations. The environment is activated at the beginning of the script and deactivated at the end.

## Data and Directory Structure

The script assumes a specific directory structure for input data and output files.  The main directory is defined by the `main_dir` variable, and the data directory is defined by the `data` variable. The script creates several subdirectories within the `data` directory to organize the results of each step.  It expects the raw VCF file (`raw_taxa.vcf.gz`) and rename chromosome text file (`rename_chrom.txt`) to be located in the `data/01_raw_data` directory.
The initial file structure looks like this.
```bash
.
├── bin
├── data
│   ├── 01_raw_data
│   │   ├── raw_taxa.vcf
│   │   ├── raw_taxa.vcf.gz
│   │   └── rename_chrom.txt
│   ├── 02_processed_data
│   │   └── subset_filter
│   ├── 03_plink
│   │   ├── not_pruned
│   │   └── pruned
│   ├── 04_pca
│   ├── 05_admixture
│   │   └── PQ_files
│   ├── 06_fst
│   └── 07_phylogeny
├── environment.yml
├── README.md
└── scripts
    ├── plot_admixture.R
    ├── plot_dendrogram.R
    ├── plot_phylo_bootstrap.R
    ├── run_admixture_parallel.py
    ├── variant_based_statistics.R
    ├── variant_processing.sh
    └── vcf2phylip.py
```


## Analysis Workflow

The `variant_processing.sh` script performs the following steps:

1. **VCF File Modification:**
    - Renames chromosome "chrZ" to "chr101" using `bcftools annotate` to ensure compatibility with PLINK and ADMIXTURE.
    - Removes the outgroup sample ("Naxos2") using `bcftools view` to avoid potential biases in downstream analyses.

```bash
#!/bin/bash

# Script to perfrom downstream analysis for VCF file

data="/home/inf-21-2024/binp28/project_binp28/version_2/data"
main_dir="/home/inf-21-2024/binp28/project_binp28/version_2"
#results=""

# Activate the conda environmnt
source activate vcf_env
echo "Activating Environment"

# Since we will be using PLINK and ADMIXTURE for our analysis. The software is not compatible with chromosomes bring letters.
# Therefore, we will change the chromosome Z (chrZ) to an aribitary chromosome 101.
# Step1: Modify the vcf.gz file
echo "Renaming chrZ to chr101"
bcftools annotate --rename-chrs ${data}/01_raw_data/rename_chrom.txt -o ${data}/01_raw_data/fixed_raw_taxa.vcf.gz \
    -O z ${data}/01_raw_data/raw_taxa.vcf.gz 

# Step2: Remove the outgroup
echo "Removing the outgroup"
bcftools view -o ${data}/01_raw_data/fixed_raw_taxa_no_outgroup.vcf.gz -O z -s ^Naxos2 ${data}/01_raw_data/fixed_raw_taxa.vcf.gz
```

2. **Variant Separation:**
    - Extracts SNPs (Single Nucleotide Polymorphisms) and INDELs (Insertions/Deletions) into separate VCF files using `bcftools view`.

```bash
# Step3: Separate the INDELS and SNPS
echo "Extracting the SNPs and INDELS in the variant calling file (VCF)."
# Extract SNPs
bcftools view -v snps -o ${data}/02_processed_data/snps_only.vcf.gz -O z ${data}/01_raw_data/fixed_raw_taxa_no_outgroup.vcf.gz 
# Extract Indels
bcftools view -v indels -o ${data}/02_processed_data/indels_only.vcf.gz -O z ${data}/01_raw_data/fixed_raw_taxa_no_outgroup.vcf.gz
```


3. **SNP Filtering:**
    - Removes multi-allelic SNPs using `bcftools view`.
    - Filters SNPs based on missing data using `vcftools`.  SNPs with a missingness rate greater than 0.999 are removed.
    - Further filters SNPs to retain only those where *no* samples have missing data.
```bash
# Step4: Remove the multi-allelic SNPs and missing SNPs from the data
echo "Removing multi-allelic SNPs"
bcftools view -m2 -M2 -v snps -o ${data}/02_processed_data/biallelic_SNPs.vcf.gz -O z ${data}/02_processed_data/snps_only.vcf.gz

echo "Filtering out missing SNPs"
vcftools --gzvcf ${data}/02_processed_data/biallelic_SNPs.vcf.gz --max-missing 0.999 --recode --stdout | bgzip -c > \
    ${data}/02_processed_data/biallelic_SNPs_filtered.vcf.gz

# Step5 - Retaining SNPs where the samples have no missing data
echo "Retaining SNPs where none of the samples have any missing SNP"
bcftools view ${data}/02_processed_data/biallelic_SNPs_filtered.vcf.gz --types snps -i 'F_MISSING=0' -Oz -o ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz

```
4. **Variant Subsampling and Statistics:**
    - Subsamples 1% of the biallelic SNPs using `vcfrandomsample` for initial variant statistics calculations.
    - Calculates various variant statistics using `vcftools`:
        - Allele frequency (`--freq`)
        - Mean depth per individual (`--depth`)
        - Mean depth per site (`--site-mean-depth`)
        - Site quality (`--site-quality`)
        - Proportion of missing data per individual (`--missing-indv`)
        - Proportion of missing data per site (`--missing-site`)
        - Heterozygosity and inbreeding coefficient per individual (`--het`)
    - Plots these statistics using an R script (`variant_based_statistics.R`).
```bash
# Step6 - Flitering
echo "Calculating the number of variants"
echo "The number of variants in the original raw file is $(bcftools view -H ${data}/01_raw_data/fixed_raw_taxa.vcf.gz | wc -l)"
echo "The number of variants in the biallelic SNPs that were filtered is $(bcftools view -H ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz | wc -l)"
echo "Sub-sampling from the biallelic SNPs file. Please wait..."
bcftools view ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz | vcfrandomsample -r 0.01 > ${data}/02_processed_data/subset.vcf       
bgzip ${data}/02_processed_data/subset.vcf
echo "The number of variants in the subset are $(bcftools view -H ${data}/02_processed_data/subset.vcf.gz | wc -l)"

SUBSET_VCF=${data}/02_processed_data/subset.vcf.gz
OUT=${data}/02_processed_data/subset_filter/subset_filter

# Calculate the allele frequency
echo "Calcualte the allele frequency"
vcftools --gzvcf $SUBSET_VCF --freq --out $OUT 

# Calculate mean depth per individual
echo "Calculate mean depth per individual"
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT

# Calculate mean depth per site
echo "Calculate mean depth per site"
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT

# Calculate site quality
echo "Calculate site quality for each site"
vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT

# Calculate proprtion of missing data per individual
echo "Calcualte proportion of missing data per individual"
vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT

# Calculate proportion of missing data per site
echo "Calculate proportion of missing data per site"
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

# Calculate heterozygosity and inbreeding coefficient per individual
echo "Calculate heterozygosity and inbreeding coefficient per individual"
vcftools --gzvcf $SUBSET_VCF --het --out $OUT


# Plot all the stats
echo "Plotting all the stats"
Rscript ${main_dir}/scripts/variant_based_statistics.R ${data}/02_processed_data/subset_filter/subset_filter ${data}/02_processed_data/subset_filter/stat_plots/
```

5. **Variant Filtering (Full Dataset):**
    - Applies more stringent filters to the full biallelic SNP dataset using `vcftools`:
        - Removes INDELs (`--remove-indels`)
        - Filters based on maximum missingness per SNP (`--max-missing`)
        - Filters based on minimum quality score (`--minQ`)
        - Filters based on minimum and maximum mean depth per individual (`--min-meanDP`, `--max-meanDP`)
        - Filters based on minimum and maximum depth per site (`--minDP`, `--maxDP`)
```bash
# Set filters
MISS=0.9
QUAL=30
MIN_DEPTH=6
MAX_DEPTH=25

echo "Applying filters"
vcftools --gzvcf ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz --remove-indels --max-missing $MISS --minQ $QUAL \
    --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > ${data}/02_processed_data/full_filtered.vcf.gz
bcftools index -t ${data}/02_processed_data/full_filtered.vcf.gz
echo "The number of variants in the fully filtered SNPs are $(bcftools view -H ${data}/02_processed_data/full_filtered.vcf.gz| wc -l)"
```

6. **Linkage Pruning:**
    - Performs linkage pruning using PLINK (`--indep-pairwise`) to reduce redundancy in the SNP dataset.

```bash
# Step7: Linkage pruning
echo "Performing Linkage pruning" 
plink -vcf  ${data}/02_processed_data/full_filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ${data}/02_processed_data/full_filtered
```

7. **PLINK Conversion:**
    - Converts the filtered VCF file to PLINK format (`.bed`, `.bim`, `.fam`) with and without incorporating the pruned data using `plink`.  The chromosome name "chr101" is changed back to "101" in the `.bim` file.

```bash
# Step8: Convert the snps to plink format by incorporating the pruned data
echo "Convert the vcf file to plink format by incorporating linkage pruning"
plink --vcf ${data}/02_processed_data/full_filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract ${data}/02_processed_data/full_filtered.prune.in --geno 0.999 --make-bed --out ${data}/03_plink/pruned/dataset_pruned
sed -i 's/^chr101/101/' ${data}/03_plink/pruned/dataset_pruned.bim


# Step9: Convert the snps file to plink format without incorporating the pruned data
echo "Convert the vcf file to plink format without incorporating linkage pruning"
plink --vcf ${data}/02_processed_data/full_filtered.vcf.gz --geno 0.999 --make-bed --out ${data}/03_plink/not_pruned/dataset --allow-extra-chr 
sed -i 's/^chr101/101/' ${data}/03_plink/not_pruned/dataset.bim
```

8. **ADMIXTURE Analysis:**
    - Runs ADMIXTURE for K values ranging from 2 to 7 using a Python script (`run_admixture_parallel.py`) that leverages the Ray framework for parallel processing. This step infers population structure.
    - Organizes the ADMIXTURE output files (`.P` and `.Q` files).
    - Collects cross-validation errors from the ADMIXTURE log files.
    - Extracts species names from the PLINK `.nosex` file for plotting.
    - Generates admixture plots using an R script (`plot_admixture.R`).
```bash
# Step10: Run ADMIXTURE
# Here try introduced parallel processing using ray framework
echo "Running ADMIXTURE on the plink file"
python3 ${main_dir}/scripts/run_admixture_parallel.py --bed_file ${data}/03_plink/pruned/dataset_pruned.bed --log_dir ${data}/05_admixture/log${i}.out --min_k 2 --max_k 7

# Move all the .P and .Q files generated by admixture analysis to the admiture directory
echo "Reorganizing files"
mv ${main_dir}/scripts/dataset_pruned* ${data}/05_admixture/PQ_files

# Get all the Cross-validation errors
echo "Get all the cross-validation erros for admixture analysis"
cat ${data}/05_admixture/log.out/*out | grep "CV" | awk '{print $3,$4}'| cut -c 4,7-20 > ${data}/05_admixture/dataset_pruned.cv.error

# Extract the species name from the individual names from the dataset_filter.nosex file.
echo "Extracting the species name for plotting"
awk '{print $2}' ${data}/03_plink/pruned/dataset_pruned.nosex > ${data}/03_plink/pruned/dataset_pruned.list


# Plotting the admixture analysis
echo "Plotting admixuture plot"
Rscript ${main_dir}/scripts/plot_admixture.R ${data}/05_admixture/PQ_files/dataset_pruned ${data}/03_plink/pruned/dataset_pruned.list 2 8 ${data}/05_admixture/admixture_plot.png 2

```

9. **Phylogenetic Analysis:**
    - Adds the outgroup sample back to the filtered VCF file for phylogenetic analysis.  Duplicate samples are removed to prevent issues with merging.
    - Converts the VCF file with the outgroup to PLINK format.
    - Runs phylogenetic analysis using two methods:
        - **Maximum Likelihood:** Uses `iqtree2` to construct a maximum likelihood tree with bootstrap support.
        - **Neighbor-Joining:** Uses `plink` to calculate pairwise distances and `fastme` to build a neighbor-joining tree.
    - Plots the phylogenetic trees using R scripts (`plot_dendrogram.R` and `plot_phylo_bootstrap.R`).
```bash
# Step 11: Phylogenetic analysis

echo "Adding the outgroup back"
# Extract chromosome and position information from the filtered VCF
bcftools query -f '%CHROM\t%POS\n' ${data}/02_processed_data/full_filtered.vcf.gz > ${data}/02_processed_data/full_filtered_positions.txt

# Use the extracted positions to filter the original VCF and retain only relevant variants for the outgroup
bcftools view -T ${data}/02_processed_data/full_filtered_positions.txt -o ${data}/02_processed_data/outgroup_filtered.vcf.gz -O z ${data}/01_raw_data/fixed_raw_taxa.vcf.gz

# Index the new outgroup VCF
bcftools index -t ${data}/02_processed_data/outgroup_filtered.vcf.gz

# Extract sample names from the filtered dataset and the outgroup dataset
bcftools query -l ${data}/02_processed_data/full_filtered.vcf.gz > ${data}/02_processed_data/full_filtered_samples.txt
bcftools query -l ${data}/02_processed_data/outgroup_filtered.vcf.gz > ${data}/02_processed_data/outgroup_samples.txt

# Identify duplicate sample names (i.e., those present in both datasets)
comm -12 <(sort ${data}/02_processed_data/full_filtered_samples.txt) <(sort ${data}/02_processed_data/outgroup_samples.txt) > ${data}/02_processed_data/duplicate_samples.txt

# Remove duplicate samples from the outgroup dataset
bcftools view -s ^$(paste -sd, ${data}/02_processed_data/duplicate_samples.txt) -o ${data}/02_processed_data/outgroup_filtered_no_duplicates.vcf.gz -O z ${data}/02_processed_data/outgroup_filtered.vcf.gz

# Index the new outgroup VCF (with duplicates removed)
bcftools index -t ${data}/02_processed_data/outgroup_filtered_no_duplicates.vcf.gz

# Merge the main dataset with the filtered outgroup dataset
bcftools merge -o ${data}/02_processed_data/full_filtered_with_outgroup.vcf.gz -O z \
    ${data}/02_processed_data/full_filtered.vcf.gz \
    ${data}/02_processed_data/outgroup_filtered_no_duplicates.vcf.gz


echo "Convert the vcf file to plink format without incorporating linkage pruning for phylogenetic analysis"
# Convert VCF to PLINK format for phylogenetic analysis
plink --vcf ${data}/02_processed_data/full_filtered_with_outgroup.vcf.gz --geno 0.999 --make-bed --out ${data}/03_plink/dataset_phylogeny --allow-extra-chr 

# Adjust chromosome name format in the .bim file (specific to dataset requirements)
sed -i 's/^chr101/101/' ${data}/03_plink/dataset_phylogeny.bim


echo "Running phylogenetic analysis by Maximum Likelihood"
# Convert VCF to Phylip format for phylogenetic analysis
python3 ${main_dir}/scripts/vcf2phylip.py -i ${data}/02_processed_data/full_filtered_with_outgroup.vcf.gz --fasta --output-folder ${data}/07_phylogeny

# Run Maximum Likelihood phylogenetic tree inference with IQ-TREE2
iqtree2 -s  ${data}/07_phylogeny/full_filtered_with_outgroup.min4.phy -m LG+I -bb 1000 -alrt 1000 -nt 64 -pre ${data}/07_phylogeny/aligned_sequences_phylogeny


echo "Running Phylogenetic analysis by Neighbour joining method"
# Compute pairwise genetic distances using PLINK
plink --bfile ${data}/03_plink/dataset_phylogeny --distance square --allow-extra-chr --out ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist

# Format sample IDs for tree construction
awk -F'\t' '{print $1"_"$2}' ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist.id > ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist.id_joined
paste -d "\t" <(cut -f1 ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist.id_joined) ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist > ${data}/07_phylogeny/phylo.dist

# Format distance matrix for FastME tree reconstruction
echo "$(wc -l < ${data}/07_phylogeny/phylo.dist)" > ${data}/07_phylogeny/phylo_formatted.dist && cat ${data}/07_phylogeny/phylo.dist >> ${data}/07_phylogeny/phylo_formatted.dist
awk '{ for(i=2;i<=NF;i++) $i=sprintf("%.6f", $i); print }' ${data}/07_phylogeny/phylo_formatted.dist > ${data}/07_phylogeny/phylo_formatted_fixed.dist

# Construct a phylogenetic tree using FastME
fastme -i ${data}/07_phylogeny/phylo_formatted_fixed.dist -o ${data}/07_phylogeny/phylo_tree_output.nwk

# Format tree output to remove redundant duplicate names
sed -E 's/\b([A-Za-z0-9]+)_\1\b/\1/g' ${data}/07_phylogeny/phylo_tree_output.nwk > ${data}/07_phylogeny/phylo_tree_formatted_output.nwk


echo "Plotting phylogenetic tree without bootstrapping"
# Plot phylogenetic tree using R script (without bootstrap support)
Rscript ${main_dir}/scripts/plot_dendrogram.R ${data}/07_phylogeny/phylo_tree_formatted_output.nwk "Naxos2" ${data}/07_phylogeny/phylogeny.png


echo "Plotting phylogenetic tree with bootstrapping"
# Plot bootstrapped phylogenetic tree using R script
Rscript ${main_dir}/scripts/plot_phylo_bootstrap.R ${data}/07_phylogeny/aligned_sequences_phylogeny.treefile "Naxos2" ${data}/07_phylogeny/phylogeny_bootstrap.png

```

10. **Environment Deactivation:**
    - Deactivates the conda environment.
```bash
# Deactivate environment
echo "Deactivating Environment"
conda deactivate
echo "Done!. Finished all analysis"
```

## Reproducibility

The entire analysis is encapsulated in a single bash script, making the workflow highly reproducible.  By running the script with the appropriate input data, the same results can be obtained.

```bash
#!/bin/bash

# Script to perfrom downstream analysis for VCF file

data="/home/inf-21-2024/binp28/project_binp28/version_2/data"
main_dir="/home/inf-21-2024/binp28/project_binp28/version_2"
#results=""

# Activate the conda environmnt
source activate vcf_env
echo "Activating Environment"

# Since we will be using PLINK and ADMIXTURE for our analysis. The software is not compatible with chromosomes bring letters.
# Therefore, we will change the chromosome Z (chrZ) to an aribitary chromosome 101.

# Step1: Modify the vcf.gz file
echo "Renaming chrZ to chr101"
bcftools annotate --rename-chrs ${data}/01_raw_data/rename_chrom.txt -o ${data}/01_raw_data/fixed_raw_taxa.vcf.gz \
    -O z ${data}/01_raw_data/raw_taxa.vcf.gz 

# Step2: Remove the outgroup
echo "Removing the outgroup"
bcftools view -o ${data}/01_raw_data/fixed_raw_taxa_no_outgroup.vcf.gz -O z -s ^Naxos2 ${data}/01_raw_data/fixed_raw_taxa.vcf.gz


# Step3: Separate the INDELS and SNPS
echo "Extracting the SNPs and INDELS in the variant calling file (VCF)."
# Extract SNPs
bcftools view -v snps -o ${data}/02_processed_data/snps_only.vcf.gz -O z ${data}/01_raw_data/fixed_raw_taxa_no_outgroup.vcf.gz 
# Extract Indels
bcftools view -v indels -o ${data}/02_processed_data/indels_only.vcf.gz -O z ${data}/01_raw_data/fixed_raw_taxa_no_outgroup.vcf.gz


# Step4: Remove the multi-allelic SNPs and missing SNPs from the data
echo "Removing multi-allelic SNPs"
bcftools view -m2 -M2 -v snps -o ${data}/02_processed_data/biallelic_SNPs.vcf.gz -O z ${data}/02_processed_data/snps_only.vcf.gz

echo "Filtering out missing SNPs"
vcftools --gzvcf ${data}/02_processed_data/biallelic_SNPs.vcf.gz --max-missing 0.999 --recode --stdout | bgzip -c > \
    ${data}/02_processed_data/biallelic_SNPs_filtered.vcf.gz

# Step5 - Retaining SNPs where the samples have no missing data
echo "Retaining SNPs where none of the samples have any missing SNP"
bcftools view ${data}/02_processed_data/biallelic_SNPs_filtered.vcf.gz --types snps -i 'F_MISSING=0' -Oz -o ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz


# Step6 - Flitering
echo "Calculating the number of variants"
echo "The number of variants in the original raw file is $(bcftools view -H ${data}/01_raw_data/fixed_raw_taxa.vcf.gz | wc -l)"
echo "The number of variants in the biallelic SNPs that were filtered is $(bcftools view -H ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz | wc -l)"
echo "Sub-sampling from the biallelic SNPs file. Please wait..."
bcftools view ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz | vcfrandomsample -r 0.01 > ${data}/02_processed_data/subset.vcf       
bgzip ${data}/02_processed_data/subset.vcf
echo "The number of variants in the subset are $(bcftools view -H ${data}/02_processed_data/subset.vcf.gz | wc -l)"

SUBSET_VCF=${data}/02_processed_data/subset.vcf.gz
OUT=${data}/02_processed_data/subset_filter/subset_filter

# Calculate the allele frequency
echo "Calcualte the allele frequency"
vcftools --gzvcf $SUBSET_VCF --freq --out $OUT 

# Calculate mean depth per individual
echo "Calculate mean depth per individual"
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT

# Calculate mean depth per site
echo "Calculate mean depth per site"
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT

# Calculate site quality
echo "Calculate site quality for each site"
vcftools --gzvcf $SUBSET_VCF --site-quality --out $OUT

# Calculate proprtion of missing data per individual
echo "Calcualte proportion of missing data per individual"
vcftools --gzvcf $SUBSET_VCF --missing-indv --out $OUT

# Calculate proportion of missing data per site
echo "Calculate proportion of missing data per site"
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT

# Calculate heterozygosity and inbreeding coefficient per individual
echo "Calculate heterozygosity and inbreeding coefficient per individual"
vcftools --gzvcf $SUBSET_VCF --het --out $OUT


# Plot all the stats
echo "Plotting all the stats"
Rscript ${main_dir}/scripts/variant_based_statistics.R ${data}/02_processed_data/subset_filter/subset_filter ${data}/02_processed_data/subset_filter/stat_plots/


# Set filters
MISS=0.9
QUAL=30
MIN_DEPTH=6
MAX_DEPTH=25

echo "Applying filters"
vcftools --gzvcf ${data}/02_processed_data/biallelic_SNPs_filtered_pre_filtering.vcf.gz --remove-indels --max-missing $MISS --minQ $QUAL \
    --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | bgzip -c > ${data}/02_processed_data/full_filtered.vcf.gz
bcftools index -t ${data}/02_processed_data/full_filtered.vcf.gz
echo "The number of variants in the fully filtered SNPs are $(bcftools view -H ${data}/02_processed_data/full_filtered.vcf.gz| wc -l)"


# Step7: Linkage pruning
echo "Performing Linkage pruning" 
plink -vcf  ${data}/02_processed_data/full_filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ${data}/02_processed_data/full_filtered


# Step8: Convert the snps to plink format by incorporating the pruned data
echo "Convert the vcf file to plink format by incorporating linkage pruning"
plink --vcf ${data}/02_processed_data/full_filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract ${data}/02_processed_data/full_filtered.prune.in --geno 0.999 --make-bed --out ${data}/03_plink/pruned/dataset_pruned
sed -i 's/^chr101/101/' ${data}/03_plink/pruned/dataset_pruned.bim


# Step9: Convert the snps file to plink format without incorporating the pruned data
echo "Convert the vcf file to plink format without incorporating linkage pruning"
plink --vcf ${data}/02_processed_data/full_filtered.vcf.gz --geno 0.999 --make-bed --out ${data}/03_plink/not_pruned/dataset --allow-extra-chr 
sed -i 's/^chr101/101/' ${data}/03_plink/not_pruned/dataset.bim


# Step10: Run ADMIXTURE

# Here try introduced parallel processing using ray framework
echo "Running ADMIXTURE on the plink file"
python3 ${main_dir}/scripts/run_admixture_parallel.py --bed_file ${data}/03_plink/pruned/dataset_pruned.bed --log_dir ${data}/05_admixture/log${i}.out --min_k 2 --max_k 7

# Move all the .P and .Q files generated by admixture analysis to the admiture directory
echo "Reorganizing files"
mv ${main_dir}/scripts/dataset_pruned* ${data}/05_admixture/PQ_files

# Get all the Cross-validation errors
echo "Get all the cross-validation erros for admixture analysis"
cat ${data}/05_admixture/log.out/*out | grep "CV" | awk '{print $3,$4}'| cut -c 4,7-20 > ${data}/05_admixture/dataset_pruned.cv.error

# Extract the species name from the individual names from the dataset_filter.nosex file.
echo "Extracting the species name for plotting"
awk '{print $2}' ${data}/03_plink/pruned/dataset_pruned.nosex > ${data}/03_plink/pruned/dataset_pruned.list


# Plotting the admixture analysis
echo "Plotting admixuture plot"
Rscript ${main_dir}/scripts/plot_admixture.R ${data}/05_admixture/PQ_files/dataset_pruned ${data}/03_plink/pruned/dataset_pruned.list 2 8 ${data}/05_admixture/admixture_plot.png 2


# Step 11: Phylogenetic analysis

echo "Adding the outgroup back"
# Extract chromosome and position information from the filtered VCF
bcftools query -f '%CHROM\t%POS\n' ${data}/02_processed_data/full_filtered.vcf.gz > ${data}/02_processed_data/full_filtered_positions.txt

# Use the extracted positions to filter the original VCF and retain only relevant variants for the outgroup
bcftools view -T ${data}/02_processed_data/full_filtered_positions.txt -o ${data}/02_processed_data/outgroup_filtered.vcf.gz -O z ${data}/01_raw_data/fixed_raw_taxa.vcf.gz

# Index the new outgroup VCF
bcftools index -t ${data}/02_processed_data/outgroup_filtered.vcf.gz

# Extract sample names from the filtered dataset and the outgroup dataset
bcftools query -l ${data}/02_processed_data/full_filtered.vcf.gz > ${data}/02_processed_data/full_filtered_samples.txt
bcftools query -l ${data}/02_processed_data/outgroup_filtered.vcf.gz > ${data}/02_processed_data/outgroup_samples.txt

# Identify duplicate sample names (i.e., those present in both datasets)
comm -12 <(sort ${data}/02_processed_data/full_filtered_samples.txt) <(sort ${data}/02_processed_data/outgroup_samples.txt) > ${data}/02_processed_data/duplicate_samples.txt

# Remove duplicate samples from the outgroup dataset
bcftools view -s ^$(paste -sd, ${data}/02_processed_data/duplicate_samples.txt) -o ${data}/02_processed_data/outgroup_filtered_no_duplicates.vcf.gz -O z ${data}/02_processed_data/outgroup_filtered.vcf.gz

# Index the new outgroup VCF (with duplicates removed)
bcftools index -t ${data}/02_processed_data/outgroup_filtered_no_duplicates.vcf.gz

# Merge the main dataset with the filtered outgroup dataset
bcftools merge -o ${data}/02_processed_data/full_filtered_with_outgroup.vcf.gz -O z \
    ${data}/02_processed_data/full_filtered.vcf.gz \
    ${data}/02_processed_data/outgroup_filtered_no_duplicates.vcf.gz


echo "Convert the vcf file to plink format without incorporating linkage pruning for phylogenetic analysis"
# Convert VCF to PLINK format for phylogenetic analysis
plink --vcf ${data}/02_processed_data/full_filtered_with_outgroup.vcf.gz --geno 0.999 --make-bed --out ${data}/03_plink/dataset_phylogeny --allow-extra-chr 

# Adjust chromosome name format in the .bim file (specific to dataset requirements)
sed -i 's/^chr101/101/' ${data}/03_plink/dataset_phylogeny.bim


echo "Running phylogenetic analysis by Maximum Likelihood"
# Convert VCF to Phylip format for phylogenetic analysis
python3 ${main_dir}/scripts/vcf2phylip.py -i ${data}/02_processed_data/full_filtered_with_outgroup.vcf.gz --fasta --output-folder ${data}/07_phylogeny

# Run Maximum Likelihood phylogenetic tree inference with IQ-TREE2
iqtree2 -s  ${data}/07_phylogeny/full_filtered_with_outgroup.min4.phy -m LG+I -bb 1000 -alrt 1000 -nt 64 -pre ${data}/07_phylogeny/aligned_sequences_phylogeny


echo "Running Phylogenetic analysis by Neighbour joining method"
# Compute pairwise genetic distances using PLINK
plink --bfile ${data}/03_plink/dataset_phylogeny --distance square --allow-extra-chr --out ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist

# Format sample IDs for tree construction
awk -F'\t' '{print $1"_"$2}' ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist.id > ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist.id_joined
paste -d "\t" <(cut -f1 ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist.id_joined) ${data}/07_phylogeny/dataset_phylogeny_SNPs_dist.dist > ${data}/07_phylogeny/phylo.dist

# Format distance matrix for FastME tree reconstruction
echo "$(wc -l < ${data}/07_phylogeny/phylo.dist)" > ${data}/07_phylogeny/phylo_formatted.dist && cat ${data}/07_phylogeny/phylo.dist >> ${data}/07_phylogeny/phylo_formatted.dist
awk '{ for(i=2;i<=NF;i++) $i=sprintf("%.6f", $i); print }' ${data}/07_phylogeny/phylo_formatted.dist > ${data}/07_phylogeny/phylo_formatted_fixed.dist

# Construct a phylogenetic tree using FastME
fastme -i ${data}/07_phylogeny/phylo_formatted_fixed.dist -o ${data}/07_phylogeny/phylo_tree_output.nwk

# Format tree output to remove redundant duplicate names
sed -E 's/\b([A-Za-z0-9]+)_\1\b/\1/g' ${data}/07_phylogeny/phylo_tree_output.nwk > ${data}/07_phylogeny/phylo_tree_formatted_output.nwk


echo "Plotting phylogenetic tree without bootstrapping"
# Plot phylogenetic tree using R script (without bootstrap support)
Rscript ${main_dir}/scripts/plot_dendrogram.R ${data}/07_phylogeny/phylo_tree_formatted_output.nwk "Naxos2" ${data}/07_phylogeny/phylogeny.png


echo "Plotting phylogenetic tree with bootstrapping"
# Plot bootstrapped phylogenetic tree using R script
Rscript ${main_dir}/scripts/plot_phylo_bootstrap.R ${data}/07_phylogeny/aligned_sequences_phylogeny.treefile "Naxos2" ${data}/07_phylogeny/phylogeny_bootstrap.png


# Deactivate environment
echo "Deactivating Environment"
conda deactivate

echo "Done!. Finished all analysis"
```

## Usage

To run the script, you will need to:

1. Install the required software packages and create the `vcf_env` conda environment.
```bash
conda env create -f environment.yml
```
2. Place your raw VCF file (`raw_taxa.vcf.gz`) and rename chromosome text file (`rename_chrom.txt`) in the `data/01_raw_data` directory. Also create the initial file structure as it is. This ensures that the scripts runs without any problem.
3. Modify the script variables (`data`, `main_dir`) if necessary to match your directory structure.
4. Execute the script: `bash variant_processing.sh`
```bash
bash variant_processing.sh
```
5. You can monitor the progress of steps using the `top` or `htop` command.
6. The final output after sucessfully running the script looks like this.
```bash
.
├── bin
├── data
│   ├── 01_raw_data
│   │   ├── fixed_raw_taxa_no_outgroup.vcf.gz
│   │   ├── fixed_raw_taxa.vcf.gz
│   │   ├── raw_taxa.vcf
│   │   ├── raw_taxa.vcf.gz
│   │   └── rename_chrom.txt
│   ├── 02_processed_data
│   │   ├── biallelic_SNPs_filtered_pre_filtering.vcf.gz
│   │   ├── biallelic_SNPs_filtered.vcf.gz
│   │   ├── biallelic_SNPs.vcf.gz
│   │   ├── duplicate_samples.txt
│   │   ├── full_filtered.log
│   │   ├── full_filtered.nosex
│   │   ├── full_filtered_positions.txt
│   │   ├── full_filtered.prune.in
│   │   ├── full_filtered.prune.out
│   │   ├── full_filtered_samples.txt
│   │   ├── full_filtered.vcf.gz
│   │   ├── full_filtered.vcf.gz.tbi
│   │   ├── full_filtered_with_outgroup.vcf.gz
│   │   ├── indels_only.vcf.gz
│   │   ├── outgroup_filtered_no_duplicates.vcf.gz
│   │   ├── outgroup_filtered_no_duplicates.vcf.gz.tbi
│   │   ├── outgroup_filtered.vcf.gz
│   │   ├── outgroup_filtered.vcf.gz.tbi
│   │   ├── outgroup_samples.txt
│   │   ├── snps_only.vcf.gz
│   │   ├── subset_filter
│   │   │   ├── stat_plots
│   │   │   ├── subset_filter.frq
│   │   │   ├── subset_filter.het
│   │   │   ├── subset_filter.idepth
│   │   │   ├── subset_filter.imiss
│   │   │   ├── subset_filter.ldepth.mean
│   │   │   ├── subset_filter.lmiss
│   │   │   ├── subset_filter.log
│   │   │   └── subset_filter.lqual
│   │   └── subset.vcf.gz
│   ├── 03_plink
│   │   ├── dataset_phylogeny.bed
│   │   ├── dataset_phylogeny.bim
│   │   ├── dataset_phylogeny.fam
│   │   ├── dataset_phylogeny.log
│   │   ├── dataset_phylogeny.nosex
│   │   ├── not_pruned
│   │   │   ├── dataset.bed
│   │   │   ├── dataset.bim
│   │   │   ├── dataset.fam
│   │   │   ├── dataset.log
│   │   │   └── dataset.nosex
│   │   └── pruned
│   │       ├── dataset_pruned.bed
│   │       ├── dataset_pruned.bim
│   │       ├── dataset_pruned.fam
│   │       ├── dataset_pruned.list
│   │       ├── dataset_pruned.log
│   │       └── dataset_pruned.nosex
│   ├── 04_pca
│   ├── 05_admixture
│   │   ├── admixture_plot_individual_k_2.png
│   │   ├── admixture_plot.png
│   │   ├── dataset_pruned.cv.error
│   │   ├── log.out
│   │   │   ├── log2.out
│   │   │   ├── log3.out
│   │   │   ├── log4.out
│   │   │   ├── log5.out
│   │   │   ├── log6.out
│   │   │   └── log7.out
│   │   └── PQ_files
│   │       ├── dataset_pruned.2.P
│   │       ├── dataset_pruned.2.Q
│   │       ├── dataset_pruned.3.P
│   │       ├── dataset_pruned.3.Q
│   │       ├── dataset_pruned.4.P
│   │       ├── dataset_pruned.4.Q
│   │       ├── dataset_pruned.5.P
│   │       ├── dataset_pruned.5.Q
│   │       ├── dataset_pruned.6.P
│   │       ├── dataset_pruned.6.Q
│   │       ├── dataset_pruned.7.P
│   │       └── dataset_pruned.7.Q
│   ├── 06_fst
│   └── 07_phylogeny
│       ├── aligned_sequences_phylogeny.bionj
│       ├── aligned_sequences_phylogeny.ckp.gz
│       ├── aligned_sequences_phylogeny.contree
│       ├── aligned_sequences_phylogeny.iqtree
│       ├── aligned_sequences_phylogeny.log
│       ├── aligned_sequences_phylogeny.mldist
│       ├── aligned_sequences_phylogeny.splits.nex
│       ├── aligned_sequences_phylogeny.treefile
│       ├── dataset_phylogeny_SNPs_dist.dist
│       ├── dataset_phylogeny_SNPs_dist.dist.id
│       ├── dataset_phylogeny_SNPs_dist.dist.id_joined
│       ├── dataset_phylogeny_SNPs_dist.log
│       ├── dataset_phylogeny_SNPs_dist.nosex
│       ├── full_filtered_with_outgroup.min4.fasta
│       ├── full_filtered_with_outgroup.min4.phy
│       ├── phylo.dist
│       ├── phylo_formatted.dist
│       ├── phylo_formatted_fixed.dist
│       ├── phylo_formatted_fixed.dist_fastme_stat.txt
│       ├── phylogeny_bootstrap.png
│       ├── phylogeny.png
│       ├── phylo_tree_formatted_output.nwk
│       └── phylo_tree_output.nwk
├── environment.yml
├── README.md
└── scripts
    ├── out.log
    ├── plot_admixture.R
    ├── plot_dendrogram.R
    ├── plot_phylo_bootstrap.R
    ├── run_admixture_parallel.py
    ├── variant_based_statistics.R
    ├── variant_processing.sh
    └── vcf2phylip.py

```

## Helper scripts

1. Prallel processing of ADMIXTURE analysis using ray framework.
```python
import subprocess
import os
import argparse
import ray


def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ADMIXTURE in parallel using Ray.")
    parser.add_argument("--bed_file", required=True, help="Path to the PLINK BED file.")
    parser.add_argument("--log_dir", required=True, help="Directory to store log files.")
    parser.add_argument("--min_k", type=int, required=True, help="Minimum K value.")
    parser.add_argument("--max_k", type=int, required=True, help="Maximum K value.")
    return parser.parse_args()

# Initialize Ray with all available resources
ray.init(num_cpus=64) #num_cpus=os.cpu_count()

@ray.remote
def run_admixture(k, bed_file, log_dir):
    #Runs ADMIXTURE for a given k in parallel.
    log_file = os.path.join(log_dir, f"log{k}.out")
    cmd = f"admixture --cv {bed_file} {k} > {log_file}"
    subprocess.run(cmd, shell=True, check=True)
    return f"Finished ADMIXTURE for K={k}"

if __name__ == "__main__":
    args = parse_arguments()
    
    # Ensure log directory exists
    os.makedirs(args.log_dir, exist_ok=True)
    
    k_values = list(range(args.min_k, args.max_k + 1))
    futures = [run_admixture.remote(k, args.bed_file, args.log_dir) for k in k_values]
    
    # Wait for all tasks to complete
    results = ray.get(futures)
    
    for res in results:
        print(res)
    
    # Shutdown Ray
    ray.shutdown()
```
2. Script to plot the admixture ananlysis
```bash
#!/usr/bin/env Rscript

rm(list=ls())

# List of necessary packages
packages <- c("ggplot2", "tidyverse", "dplyr")

# Install packages if they are not already installed
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript admixture_plot.R <file_prefix> <sample_path> <min_k> <max_k> <img_path> <individual_k>")
}

file_prefix = args[1]  # Example: "data/05_admixture/PQ_files/dataset_filter"
samplelist_path = args[2]
min_k = as.numeric(args[3])  # Example: 3
max_k = as.numeric(args[4])  # Example: 5
save_img = args[5]
individual_k = as.numeric(args[6])  # The specific k-value for individual plot

# Check for valid input
if (is.na(min_k) | is.na(max_k) | min_k > max_k) {
  stop("Invalid K range. Ensure min_k and max_k are numeric and min_k <= max_k.")
}

if (is.na(individual_k) | individual_k < min_k | individual_k > max_k) {
  stop("Invalid individual_k value. Ensure individual_k is between min_k and max_k.")
}

# Read sample list
if (!file.exists(samplelist_path)) {
  stop(paste("Sample list file not found:", samplelist_path))
}

samplelist = read_tsv(samplelist_path, col_names = "sample", show_col_types = FALSE)

# Initialize tibble to store all data
all_data = tibble(sample=character(), k=numeric(), Q=character(), value=numeric())

# Loop over K values to read data
for (k in seq(min_k, max_k)) {
  file_path = paste0(file_prefix, ".", k, ".Q")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found, skipping:", file_path))
    next
  }
  
  # Read Q matrix
  data = read_delim(file_path, 
                    col_names = paste0("Q", seq(1:k)), 
                    delim = " ", 
                    show_col_types = FALSE)
  
  # Add sample and k values
  if (nrow(data) != nrow(samplelist)) {
    stop(paste("Mismatch between sample list and Q matrix for k =", k))
  }
  
  data$sample = samplelist$sample
  data$k = k
  
  # Reshape to long format
  data = data %>% pivot_longer(cols = starts_with("Q"), names_to = "Q", values_to = "value")
  
  # Append to all_data
  all_data = bind_rows(all_data, data)
}

# Plot all K values with facet_wrap
p_all = ggplot(all_data, aes(x = sample, y = value, fill = factor(Q))) + 
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14)
  )+
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  facet_wrap(~k, ncol = 1)

# Save the plot with all K values
ggsave(save_img, plot = p_all, width = 12, height = 6, dpi = 300, device = "png")


# Extract directory and base name from save_img
save_dir <- dirname(save_img)
save_base <- tools::file_path_sans_ext(basename(save_img))

# Plot for individual k value
p_individual = ggplot(all_data[all_data$k == individual_k, ], aes(x = sample, y = value, fill = factor(Q))) + 
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 14)
        )+
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  ggtitle(paste("Individual Plot for K =", individual_k))

# Construct the new path for the individual plot
individual_plot_path <- file.path(save_dir, paste0(save_base, "_individual_k_", individual_k, ".png"))

# Save the individual plot
ggsave(individual_plot_path, plot = p_individual, width = 12, height = 6, dpi = 300, device = "png")


print("All plots saved successfully!")
```
3. Conversion of vcf to phylip format. This is the only helper script that was borrowed from this repository Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861.

4. Script to plot dendrograms

4a. Maximum likelihood method

```bash
#!/usr/bin/env Rscript

rm(list = ls())

# Set CRAN and Bioconductor repositories
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install Bioconductor manager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List required packages
packages <- c("ape", "ggplot2", "tidyverse", "ggtree")

# Install missing packages
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "ggtree") {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

library(ape)
library(ggtree)

print("All packages loaded successfully!")

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript phylo_plot_bootstrap.R <tree_file> <outgroup> <img_path>")
}

tree_file = args[1]
outgroup = args[2]
save_img = args[3]

# Check if tree file exists
if (!file.exists(tree_file)) {
  stop(paste("Tree file not found:", tree_file))
}

# Read and root the tree (handling potential polytomies)
tree <- read.tree(tree_file)

if (!(outgroup %in% tree$tip.label)) {
  stop(paste("Outgroup", outgroup, "not found in the tree!"))
}

rooted_tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
rooted_tree <- ladderize(rooted_tree) # Improve readability



# Create the plot
png(save_img, width = 1600, height = 1000, res = 200)

p <- ggtree(rooted_tree, layout = "rectangular") +
  geom_tiplab(size = 3, hjust = -0.1) +
  theme_tree2() +
  xlim(0, max(node.depth.edgelength(rooted_tree)) * 1.2) +
  labs(title = "Rooted Phylogenetic Tree with Bootstrap Values")

# Add branch lengths (if needed)

p <- p + geom_treescale()

# Add bootstrap values to the branches (if they exist)
if (!is.null(rooted_tree$node.label) && length(rooted_tree$node.label) > 0 && !any(is.na(as.numeric(gsub("/.*", "", rooted_tree$node.label))))) {
  p <- p + geom_nodelab(aes(label = gsub("/.*", "", node.label)), size = 3, nudge_x = -0.01, nudge_y = 0.02, angle = 45)  # Adjust nudge and angle
}

p

dev.off()

print("Done!")
```

4b. Neighbour Joining method

```bash
#!/usr/bin/env Rscript

rm(list=ls())

# Set CRAN and Bioconductor repositories
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install Bioconductor manager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List required packages
packages <- c("ape", "ggplot2", "tidyverse", "ggtree")

# Install missing packages
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "ggtree") {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

library(ape)
library(ggtree)

print("All packages loaded successfully!")

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript phylo_plot.R <tree_file> <outgroup> <img_path>")
}

tree_file = args[1]
outgroup = args[2]
save_img = args[3]

# Check if tree file exists
if (!file.exists(tree_file)) {
  stop(paste("Tree file not found:", tree_file))
}

# Read and root the tree
tree <- read.tree(tree_file)
if (!(outgroup %in% tree$tip.label)) {
  stop(paste("Outgroup", outgroup, "not found in the tree!"))
}
rooted_tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
rooted_tree <- ladderize(rooted_tree)  # Improve readability

# Save the tree plot
png(save_img, width = 1600, height = 1000, res = 200)

ggtree(rooted_tree, layout = "rectangular") + 
  geom_tiplab(size = 3, hjust = -0.1) + 
  theme_tree2() +  
  xlim(0, max(node.depth.edgelength(rooted_tree)) * 1.2) +  
  labs(title = "Rooted Phylogenetic Tree")

dev.off()

print("Done!")
```

5. Varaint Based Statistics
```bash
#!/usr/bin/env Rscript

rm(list=ls())

# List of necessary packages
packages <- c("ggplot2", "tidyverse", "dplyr")

# Install packages if they are not already installed
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

library(tidyverse)
library(dplyr)
library(ggplot2)

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript variants_based_statistics.R <input_prefix> <output_directory>")
}

input_prefix = args[1]  # Example: "./data/01_raw_data/subset_filter/subset_filter"
output_dir = args[2]  # Example: "./output"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to save plot
save_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300)
}

# Read and plot quality scores
var_qual = read_delim(paste0(input_prefix, ".lqual"), delim = "\t", col_names = c("chr","pos","qual"), skip=1)
var_qual_plot = ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue1", colour="black", alpha = 0.3) +
  theme_light()
save_plot(var_qual_plot, file.path(output_dir, "subset_filter.lqual.png"))

# Read and plot depth distribution
var_depth = read_delim(paste0(input_prefix, ".ldepth.mean"), delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
var_depth_plot = ggplot(var_depth, aes(mean_depth)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(var_depth_plot, file.path(output_dir, "subset_filter.ldepth.mean.png"))

# Read and plot missingness
var_miss = read_delim(paste0(input_prefix, ".lmiss"), delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
var_miss_plot = ggplot(var_miss, aes(fmiss)) + 
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(var_miss_plot, file.path(output_dir, "subset_filter.lmiss.png"))

# Read and plot individual depth
ind_depth = read_delim(paste0(input_prefix, ".idepth"), delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)
ind_depth_plot = ggplot(ind_depth, aes(depth)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(ind_depth_plot, file.path(output_dir, "subset_filter.idepth.png"))

# Read and plot individual missingness
ind_miss = read_delim(paste0(input_prefix, ".imiss"), delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
ind_miss_plot = ggplot(ind_miss, aes(fmiss)) +
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  theme_light()
save_plot(ind_miss_plot, file.path(output_dir, "subset_filter.imiss.png"))

# Read and plot heterozygosity
ind_het = read_delim(paste0(input_prefix, ".het"), delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
ind_het_plot = ggplot(ind_het, aes(f)) + 
  geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) + 
  theme_light()
save_plot(ind_het_plot, file.path(output_dir, "subset_filter.het.png"))

print("Done!")

```
This README provides a detailed explanation of the steps involved in the downstream analysis of the VCF file, making the workflow transparent and reproducible.  The clear organization of the script and the use of a conda environment further enhance the reproducibility of the analysis.
