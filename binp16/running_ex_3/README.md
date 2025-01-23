# Genetic Analysis of Anastasia Romanov

**Authors**:  
- Jacob Daigle  
- Jyothishree Veerabhadran  
- Chandrashekar CR  

---

## Introduction

This project focuses on using genetic distance analysis to determine familial relationships within the Romanov family, specifically to investigate the true identity of Anastasia Romanov. By comparing mitochondrial DNA (mtDNA) and Y chromosome sequences, we calculate genetic distances using alignment and identity scores. These distances give insight into the genetic similarity and potential maternal inheritance patterns, offering clues toward identifying potential claimants to the Romanov lineage.

## Project Structure

The analysis is divided into three parts, with each author responsible for specific tasks:

### Part I (Student 1 - Jacob Daigle)
1. **Parser**  
   Parses a given input file containing individual names, hemophilia status, mtDNA sequences, and Y chromosome sequences. Outputs two FASTA files: one for mtDNA and one for Y chromosome sequences, each with individual names as headers.

2. **Multiple Sequence Aligner**  
   Aligns the sequences in the parsed FASTA files, calculating identity and alignment scores for each pair. The output is a tab-delimited text file with pairwise similarity scores (identity and alignment) for each pair of individuals.

### Part II (Student 2 - Jyothishree Veerabhadran)
1. **Building the Genetic Distance Matrix**  
   Reads the calculated identity and alignment scores to create genetic distance matrices. Outputs two matrices (tab-delimited): one for identity scores and one for alignment scores.

2. **Identifying Most Similar Individuals**  
   This script identifies the most genetically similar individuals by reading the genetic distance matrices. It allows users to either display all pairwise similarities or find the closest genetic match for a specific individual.

### Part III (Student 3 - Chandrashekar CR)
1. **Visualization (Dendrograms and Heatmaps)**  
   This script generates dendrograms and heatmaps for the genetic similarity matrices based on identity and alignment scores. Using hierarchical clustering, it provides a visual representation of genetic relationships. The Romanov family members are color-coded for easy identification.

## Methods

- **Alignment Score Calculation**: Unknown nucleotides (represented as `?`) are scored as gaps. When paired with known nucleotides, they are scored based on weighted averages to avoid random scoring biases.
- **Default Scoring Parameters**: 
   - Match = +3
   - Transition = -1
   - Transversion = -2
   - Gap = -1

- **Genetic Distance Calculation**:
  - **Based on Alignment Score**: 
    \[
    \text{Genetic distance} = 1 - \frac{\text{Alignment score of AB}}{\text{average(Alignment scores of AA, BB)}}
    \]
  - **Based on Identity Score**:  
    \[
    \text{Genetic distance} = 1 - \frac{\text{Identity score of AB}}{100}
    \]

## Results

### Dendrogram Analysis
The dendrograms provide insight into the genetic similarity between individuals in the Romanov family:
- **Alignment Score Dendrogram**: Shows Grigori Rasputin and Nicolai Romanov II as outgroups, clustering together. Alexandra and Tatiana Romanov also cluster closely, with Anastasia 1 and Anastasia 2 forming another distinct cluster.
- **Identity Score Dendrogram**: Shows Nicolai as an outgroup with Tatiana and Rasputin also distanced. Notably, Alexandra Romanov clusters closely with Anastasia 1, suggesting that Anastasia 1 is the likely claimant.

### Interpretation
- **mtDNA Inheritance**: mtDNA is maternally inherited, making it a powerful tool for confirming maternal lineage. Anastasia 1’s close alignment with Alexandra Romanov in the identity dendrogram suggests she may be the genuine Anastasia.
- **Scoring Strategy**: Heavier penalties on mismatches, transitions, and gaps highlight the subtle genetic differences between family members, while favoring alignments that match closely related individuals.

## Discussion

To accurately assess familial relationships in mtDNA analysis, it’s crucial to apply penalties that favor close matches and common mutations. Our scoring approach penalizes mismatches and gaps while moderately penalizing transitions, thus focusing on preserving family resemblance. This method enhances the ability to discern true genetic similarity in closely related individuals, such as family members.

## Conclusion

This project offers a comprehensive genetic analysis of the Romanov family. While genetic analysis provides valuable insights, it should be supplemented with historical records and additional evidence for accurate identification. Future studies can use similar techniques to further explore genealogical relationships among historical figures or claimants of lineage.

## Requirements

- Python 3.x
- Packages: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`

## Install Conda environment
- conda env create -f environment.yml

## File structure
- Use the file structure to undestand how the files are organized.
