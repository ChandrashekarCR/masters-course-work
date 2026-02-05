#!/bin/bash

# Script to perfrom downstream analysis for VCF file

data="/home/inf-21-2024/masters_course_work/binp28/data"
main_dir="/home/inf-21-2024/masters_course_work/binp28"
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
