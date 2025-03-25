# Malaria Phylogenetics and Comparative Genomics Pipeline

## Author
- [@ChandrashekarCR](https://github.com/ChandrashekarCR)

This repository contains a bash script (`malaria_pipeline.sh`) that processes genomic data for various *Plasmodium* species and related taxa. The workflow includes genome filtering, gene prediction, ortholog identification, phylogenetic analysis, and tree reconstruction using RAxML.

## Computational Environment
A dedicated conda environment, `malaria_env`, is used to manage the software dependencies. The script also utilizes a separate `busco_env` for BUSCO analysis. The environments are activated and deactivated as needed to ensure reproducibility and avoid conflicts.

## Data and Directory Structure
The script assumes the following directory structure for input and output data:
```bash
.
├── bin
│   ├── gffParse.pl
│   ├── gmes_petap.pl
│   └── gm_key_64
├── data
│   ├── 01_raw_data
│   ├── 02_genemark
│   ├── 03_haemoproteus
│   ├── 04_blast
│   ├── 05_proteinortho
│   ├── 06_busco
│   ├── 07_alignments
│   └── 08_raxml
├── file_struct.txt
├── README.md
├── results
└── scripts
    ├── busco_protein_fasta.py
    ├── datParser.py
    ├── filtering_scafolds.py
    └── malaria_pipeline.sh
```



## Pipeline Overview
The `malaria_pipeline.sh` script performs the following key steps:

### Step 0: Collecting Data

**Reasoning**: This step initializes the pipeline by providing the location of the required genomic data files. The data files should be stored in the `/resources/binp29/Data/malaria/` directory for further processing. 

The required genomic data files for analysis are stored at the following location on the bioinf-serv2 (Bioinformatics Course Server):

/resources/binp29/Data/malaria/

These data files include the raw genome sequence required for further processing and gene prediction tasks.

### Step 1: Data Preprocessing
1. **Computing Genome Size and GC Content**: This part calculates the genome size and GC content for each *Plasmodium* species genome file. The `grep` command is used to remove headers (`>`) in FASTA files, `tr -d "\n"` is used to concatenate sequences, and `wc -c` counts the number of bases. The GC content is calculated by counting the occurrences of `G` or `C` and dividing by the total genome size.

**Reasoning**:
- `grep -v "^>"`: Skips header lines (FASTA headers).
- `tr -d "\n"`: Removes newlines to treat the genome as a single line, so `wc -c` can count the length properly.
- `grep -o "[GCgc]" | wc -l`: Counts the number of `G` or `C` nucleotides in the genome to calculate GC content.

2. **Extracting Gene Counts**: This extracts the number of genes from the GTF files by cutting out the gene names from the 9th column and counting unique gene IDs.

**Reasoning**:
- `cut -f9`: Extracts the gene attributes from the GTF file.
- `cut -d \" -f2`: Extracts the actual gene ID value by splitting the string.
- `sort | uniq`: Sorts the gene IDs and then counts the unique ones to get the gene count.

```bash
# Step 1: Basic data analysis. Genome size, genes and genomic GC content.
# Genome size.
echo "The genome size for the plasmodium species is as follows"
for file in "$data/01_raw_data/plasmodium_data/"*.genome; do
    strain_name=$(basename $file)
    genome_size=$(cat $file | grep -v "^>" | tr -d "\n" | wc -c)
    gc_content=$(cat $file | grep -v "^>" | tr -d "\n" | grep -o "[GCgc]"| wc -l)
    gc_percentage=$(echo "($gc_content/$genome_size)*100" | bc -l)
    echo "$strain_name Genome Size: $genome_size bp, GC Content: $gc_percentage%"
done


# Genes
echo "The number of genes for the plasmodium species is as follows"
for file in "$data/02_genemark/"*.gtf; do
    strain_name=$(basename $file)
    no_of_genes=$(cat $file | cut -f9 | cut -d  \" -f2 | sort | uniq | wc -l)
    echo "$strain_name $no_of_genes"
done
```

### Step 2: Filtering Haemoproteus tartakovskyi Genome
1. **Filter Genomes Based on GC Content and Scaffold Length**: This step removes scaffolds with a GC content > 28% and length shorter than 3000 nt. The filtering script (`filtering_scafolds.py`) is used here with the specific threshold values.

**Reasoning**: 
- `--gc_content 28`: This flag sets the threshold for removing scaffolds with GC content higher than 28%.
- `--scaffold_len 30000`: This flag ensures that only scaffolds longer than 3000 nt are kept.

```bash
# Step2: Cleaning the genome sequence of Haemoproteus tartakovskyi 
# Removing all scaffolds that have a GC content greater than 28% and a length shorter than 3000 nucleotides.
echo "Removing all the scaffolds based on GC content and nucleotide content"
python3 "$scripts/filtering_scafolds.py" --fasta_file "$data/03_haemoproteus/Haemoproteus_tartakovskyi.raw.genome" --output_file "$data/03_haemoproteus/Haemoproteus_tartakovskyi_filtered.genome" --gc_content 28 --scaffold_len 30000
echo "Sequence before filtering is $(cat $data/03_haemoproteus/Haemoproteus_tartakovskyi.raw.genome | grep "^>" | wc -l)"
echo "Sequence after filtering is $(cat $data/03_haemoproteus/Haemoproteus_tartakovskyi_filtered.genome | grep "^>" | wc -l)"
```
### Step 3: Gene Prediction
1. **GeneMark Gene Prediction**: This uses the GeneMark-ES program to predict genes from the filtered genome. The `--min_contig` option ensures that only contigs of a certain minimum length are used.

**Reasoning**: 
- `--min_contig 10000`: This ensures that only contigs of at least 10,000 nt in length are considered for gene prediction, avoiding overly short or ambiguous contigs that could introduce noise.

```bash
# Step3: Gene Prediction for the new haemoproteus tartakovskyi filtered file.
mkdir $data/03_haemoproteus/haemoproteus_gene_prediction
echo "Running gene prediction on harmoproteus_tartakovskyi..."
## Min contig is 100000, is because we have the genome of the host and the parasite as well. We need to keep only the shorter contigs, but not too short ones because that would lead to no proper results.
gmes_petap.pl --ES --sequence $data/03_haemoproteus/Haemoproteus_tartakovskyi_filtered.genome --min_contig 10000 --work_dir $data/03_haemoproteus/haemoproteus_gene_prediction --cores 60
```

### Step 4: Functional Annotation with BLAST
1. **BLASTx**: The `blastx` command is used to align the predicted protein sequences against the SwissProt database to identify potential homologous genes.

**Reasoning**: 
- `-evalue 1e-5`: The E-value threshold for reporting results. This ensures that only hits with a reasonable probability of being homologous are kept.
- `-num_threads 35`: This runs BLAST using 35 threads, improving efficiency on multi-core machines.

```bash
# Formatting the gtf file for the perl script to run
cat $data/03_haemoproteus/haemoproteus_gene_prediction/haemoproteus_tartakovskyi_gene_predict.gtf | sed "s/  length=.*\tGeneMark.hmm/\tGeneMark.hmm/" > $data/03_haemoproteus/haemoproteus_gene_prediction/haemoproteus_tartakovskyi_gene_predict_formatted.gtf

# Create FASTA sequences from the gff file.
perl $bin_dir/gffParse.pl -c -p -F -i $data/03_haemoproteus/Haemoproteus_tartakovskyi_filtered.genome -g $data/03_haemoproteus/haemoproteus_gene_prediction/haemoproteus_tartakovskyi_gene_predict_formatted.gtf 
mv $scripts/gffParse* $data/03_haemoproteus/

## Step4: Blast Search
mv $data/03_haemoproteus/gffParse.fna $data/03_haemoproteus/haemoproteus_tartakovskyi_gene_predict_formatted.fna
mv $data/03_haemoproteus/gffParse.faa $data/03_haemoproteus/haemoproteus_tartakovskyi_gene_predict_formatted.faa

#echo "Running Blast.."
blastx -query $data/03_haemoproteus/haemoproteus_tartakovskyi_gene_predict_formatted.fna -db SwissProt -out $data/04_blast/haemoproteus_tartakovskyi_blast.blastx -outfmt 6 -evalue 1e-5 -num_threads 35
 Create a symlink for the taxonomy file
ln -s /resources/binp29/Data/malaria/taxonomy.dat $data/04_blast/taxonomy.dat
ln -s /resources/binp29/Data/malaria/uniprot_sprot.dat $data/04_blast/uniprot_sprot.dat

# A script is provided to retrieve the host scaffolds
echo "Script to retrieve the host scaffolds" 
python $scripts/datParser.py $data/04_blast/haemoproteus_tartakovskyi_blast.blastx $data/03_haemoproteus/haemoproteus_tartakovskyi_gene_predict_formatted.fna $data/04_blast/taxonomy.dat $data/04_blast/uniprot_sprot.dat 

```
### Step 5: Ortholog Identification
1. **Proteinortho**: The `proteinortho` command is used to identify orthologous genes across different species. It compares the protein sequences and detects orthologs that may have evolved from a common ancestor.

**Reasoning**:
- `proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa`: This command runs `proteinortho` on a predefined list of species protein sequences.

2. **BUSCO**: The BUSCO tool evaluates the completeness of a genome by searching for universal single-copy orthologs.

**Reasoning**: 
- `-m prot`: Specifies the protein mode for analysis.
- `-l apicomplexa`: Uses the specific lineage data for *Apicomplexa*.

```bash
# Step5: Proteinortho
echo "Converting the gtf files to fasta files..." 
for file in $data/01_raw_data/plasmodium_data/*; do
    strain_name=$(basename -- $file| cut -d '.' -f1 | awk '{print tolower($0)}')
    #mv $file $data/01_raw_data/plasmodium_data/$strain_name.genome
done;
cp $data/03_haemoproteus/haemoproteus_gene_prediction/haemoproteus_tartakovskyi_gene_predict_formatted.gtf $data/02_genemark/
mv $data/02_genemark/haemoproteus_tartakovskyi_gene_predict_formatted.gtf $data/02_genemark/haemoproteus_tartakovskyi.gtf
cp $data/03_haemoproteus/Haemoproteus_tartakovskyi.raw.genome $data/01_raw_data/plasmodium_data/
mv $data/01_raw_data/plasmodium_data/Haemoproteus_tartakovskyi.raw.genome $data/01_raw_data/plasmodium_data/haemoproteus_tartakovskyi.genome
mv $data/01_raw_data/plasmodium_data/plasmodium_faciparum.genome  $data/01_raw_data/plasmodium_data/plasmodium_falciparum.genome



cd $data/05_proteinortho
declare -A species_map
species_map=(
    ["plasmodium_berghei"]="Pb"
    ["plasmodium_cynomolgi"]="Pc"
    ["plasmodium_falciparum"]="Pf"
    ["plasmodium_knowlesi"]="Pk"
    ["plasmodium_vivax"]="Pv"
    ["plasmodium_yoelii"]="Py"
    ["toxoplasma_gondii"]="Tg"
    ["haemoproteus_tartakovskyi"]="Ht"
)

for file in $data/01_raw_data/plasmodium_data/*; do
    full_name=$(basename -- $file| cut -d '.' -f1 | awk '{print tolower($0)}')

    strain_name=${species_map[$full_name]:-$full_name}
    
    #mv $file $data/01_raw_data/plasmodium_data/$strain_name.genome
    #echo $data/01_raw_data/plasmodium_data/$full_name.genome
    #echo $data/02_genemark/$full_name.gtf

    perl $bin_dir/gffParse.pl -c -p -b $strain_name -F -i $data/01_raw_data/plasmodium_data/$full_name.genome -g $data/02_genemark/$full_name.gtf
    sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' $strain_name.faa
 
done

cd $scripts
cd $data/05_proteinortho
echo "Running proteinortho..."
proteinortho6.pl {Ht,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa 
cd $scripts
echo "Deactivating environment"
conda deactivate

# Step 6: BUSCO Analysis
echo "Activating BUSCO environment"
conda activate busco_env


echo "Running BUSCO analysis..."

# Define the list of species codes
species_list=(Ht Pb Pc Pf Pk Pv Py Tg)

# Loop through each species and run BUSCO
for species in "${species_list[@]}"; do

    echo "Processing $species..."
    
    # Run BUSCO
    cd $data/06_busco/
    busco -i $data/05_proteinortho/$species.faa -o $species -m prot -l apicomplexa
done
 
echo "Deactivating BUSCO environment"
conda deactivate
```
### Step 6: Multiple Sequence Alignment
1. **Clustal Omega**: Aligns protein sequences using the `clustalo` tool to generate a multiple sequence alignment.

**Reasoning**: 
- `--threads 50`: This option uses 50 threads, making the process faster on machines with many cores.
- `-o $data/07_alignments/$busconame.aligned.faa`: This specifies the output file for the aligned sequences.

```bash
# Step 7: BUSCO Filtering and Running Alignments
echo "Running BUCSO filtering..."
python3 $scripts/busco_protein_fasta.py --full_table $data/06_busco/ --fasta_directory $data/05_proteinortho/ --output_dir $data/06_busco/busco_filter_data/

echo "Running Alignments.."
for file in $data/06_busco/busco_filter_data/*.faa; do
    busconame=$(basename -- $file | cut -d "." -f1)
    clustalo -i $file -o $data/07_alignments/$busconame.aligned.faa -v --threads 50
done
```

### Step 7: Phylogenetic Tree Construction
1. **RAxML**: This step builds a phylogenetic tree using the maximum-likelihood method provided by RAxML.

**Reasoning**: 
- `-s aligned.faa`: Specifies the aligned sequence file.
- `-n tree`: Specifies the name of the output tree file.
- `-m PROTGAMMAAUTO`: Automatically selects the best substitution model for protein sequences.

```bash
# Step 8: Phylogenetic Tree Construction
echo "Running RAxML..."
for file in $data/07_alignments/*.faa; do
    aligned_name=$(basename -- "$file" | cut -d "." -f1)
    raxmlHPC -s "$file" -n "$aligned_name.tre" -o Tg -m PROTGAMMABLOSUM62 -p 12345 -w "$data/08_raxml/" -T 16
done
```

## Execution
### Create the file strucutre
Please refer the file structure and create the directories accrodingly. There might some additional folders you might create.

### Create the required enivronments using conda
```bash
conda env create -f malaria_env.yml
conda env create -f busco_env.yml
```

### Running the Pipeline
```bash
bash malaria_pipeline.sh
```

## Filtering Criteria
The following filters were applied to ensure high-quality genome sequences:
- **GC Content Threshold**: 28%
- **Minimum Scaffold Length**: 3000 bp
- **BUSCO Completeness**: Only complete or single-copy orthologs retained
- **Phylogenetic Analysis**: Single-copy orthologs used for tree reconstruction

## Results Interpretation
1. **Genome Size & GC Content**: Comparative analysis of *Plasmodium* species and related taxa.
2. **Gene Prediction & Filtering**: Removal of host scaffolds and filtering of gene predictions.
3. **Ortholog Detection**: Comparison of `proteinortho` and `BUSCO` outputs.
4. **Phylogenetic Trees**: Topology assessment and interpretation of evolutionary relationships.


## References
1. **BUSCO**: https://busco.ezlab.org/
2. **Proteinortho**: https://www.bioinf.uni-leipzig.de/Software/proteinortho/
3. **RAxML**: https://cme.h-its.org/exelixis/web/software/raxml/index.html
4. **Clustal Omega**: https://www.ebi.ac.uk/Tools/msa/clustalo/


## Note
This repository is still being updated. Any mistakes will be revised later.
