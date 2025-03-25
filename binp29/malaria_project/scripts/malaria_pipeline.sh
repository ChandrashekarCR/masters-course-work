#!/bin/bash

# This script is used for processing and filtering the plasmodium data.

data="/home/inf-21-2024/binp29/malaria_project/data"
main_dir="/home/inf-21-2024/binp29/malaria_project"
scripts="/home/inf-21-2024/binp29/malaria_project/scripts"
bin_dir="/home/inf-21-2024/binp29/malaria_project/bin"


## Source the conda.sh script to enable `conda activate` command
source ~/miniconda3/etc/profile.d/conda.sh  
## Activate your environment
echo "Activating environment"
conda activate malaria_env  # Replace `vcf_env` with your environment name
## Now you can run commands in the activated environment
echo "Environment activated"

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


# Step2: Cleaning the genome sequence of Haemoproteus tartakovskyi 
# Removing all scaffolds that have a GC content greater than 28% and a length shorter than 3000 nucleotides.
echo "Removing all the scaffolds based on GC content and nucleotide content"
python3 "$scripts/filtering_scafolds.py" --fasta_file "$data/03_haemoproteus/Haemoproteus_tartakovskyi.raw.genome" --output_file "$data/03_haemoproteus/Haemoproteus_tartakovskyi_filtered.genome" --gc_content 28 --scaffold_len 30000
echo "Sequence before filtering is $(cat $data/03_haemoproteus/Haemoproteus_tartakovskyi.raw.genome | grep "^>" | wc -l)"
echo "Sequence after filtering is $(cat $data/03_haemoproteus/Haemoproteus_tartakovskyi_filtered.genome | grep "^>" | wc -l)"


# Step3: Gene Prediction for the new haemoproteus tartakovskyi filtered file.
mkdir $data/03_haemoproteus/haemoproteus_gene_prediction
echo "Running gene prediction on harmoproteus_tartakovskyi..."
## Min contig is 100000, is because we have the genome of the host and the parasite as well. We need to keep only the shorter contigs, but not too short ones because that would lead to no proper results.
gmes_petap.pl --ES --sequence $data/03_haemoproteus/Haemoproteus_tartakovskyi_filtered.genome --min_contig 10000 --work_dir $data/03_haemoproteus/haemoproteus_gene_prediction --cores 60

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

echo "Acitvating malaria environment"
conda activate malaria_env

# Step 7: BUSCO Filtering and Running Alignments
echo "Running BUCSO filtering..."
python3 $scripts/busco_protein_fasta.py --full_table $data/06_busco/ --fasta_directory $data/05_proteinortho/ --output_dir $data/06_busco/busco_filter_data/

echo "Running Alignments.."
for file in $data/06_busco/busco_filter_data/*.faa; do
    busconame=$(basename -- $file | cut -d "." -f1)
    clustalo -i $file -o $data/07_alignments/$busconame.aligned.faa -v --threads 50
done

# Step 8: Phylogenetic Tree Construction
echo "Running RAxML..."
for file in $data/07_alignments/*.faa; do
    aligned_name=$(basename -- "$file" | cut -d "." -f1)
    raxmlHPC-PTHREADS-AVX2 -s "$file" -n "$aligned_name.tre" -o Tg -m PROTGAMMABLOSUM62 -p 12345 -w "$data/08_raxml_new/" -# 10 -T 60
done

echo "Deactivating environment"
conda deactivate
