# BINP28 Amplicon Exercise

**1. How many hindgut microbiota termite samples were studied (colony + caste)?**

The study analyzed hindgut microbiota from termites collected from eight colonies, with samples taken from different castes, including workers, soldiers, and alates (winged termites). There were 45 from colonies and 19 from caste analysis. So a total of 64 samples.

**2. What is the most abundant phylum in the termite microbiota?**

The most abundant phylum in the termite hindgut microbiota is Spirochaetes, with the genus Treponema being particularly dominant. Treponema accounted for approximately 32% of the sequences in the bacterial community

**3. How does the hindgut microbiota differ amongst castes?**

The hindgut microbiota varies significantly between castes:
- Workers and soldiers have a higher abundance of Treponema and Endomicrobia, which are crucial for lignocellulose digestion.
- Alates (winged termites) show a dramatic reduction in protist populations, such as Parabasalia and Oxymonadida, likely due to changes in diet or physiology as they prepare to swarm and establish new colonies.

**4. How much space do the files take after decompression (use disk usage du)?**

```bash
/amplicon_seq/data/AmpliconData$ du -sh */*.fastq.gz | awk '{sum+=$1} END {print sum " MB"}'

123.3 MB
```
**5. The study was done in the USA, but EBI geographical coordinates point to a different country. Which country is it? What does that tell you about the reliability of the data?**

Kyrgyzstan. The coordinates (https://www.ncbi.nlm.nih.gov/biosample/SAMEA3594527) give the location of Kyrgyzstan, but the geographical location says USA. Therefore we need to be cautious while using publicly available data.

**6. How would you remove multiple copies of a sample?**

To remove multiple copies of a sample, you can use bioinformatics tools such as:
- FastQC or MultiQC to identify duplicate sequences.
-  CD-HIT to cluster sequences and remove duplicates based on similarity
thresholds.
-  SAMtools for removing duplicate reads in sequencing data.​ These tools ensure that only unique sequences are retained for downstream analysis, improving the accuracy of the results.

**7. How many sequences do you get in total?**

```bash
for file in *.fastq; do
    cat "$file" | grep "^@M00" | wc -l
done | awk '{sum+=$1} END {print sum " sequences"}'

7850134 sequences
```

**8. What is the average read length?**

```bash
for file in *.fastq; do
    cat "$file" 
done | awk 'NR%4==2 {total+=length($0);count++} END {print "Avg length is " total/count}'

Avg length is 227.645
```

**9. What is the meaning of - - - - in the paste command?**

It splits the 4 lines of the FASTQ file.

**10. What are we doing in the second line?**

cat *.fastq | paste - - - - | cut -f 2 | grep -P "^TTAGA[CG]AC"
Here, we look for the pattern either with C or G.

**11. Study the graphs. You can compare them to this example report https://www.bioinformatics.babraham.ac.uk/projects/fastqc/? What is the meaning of two graphs?**

A comparison of the reads with the example report revealed that the sequencing run suffered from low base quality scores and significant adapter contamination, particularly towards the ends of the reads. This poor quality necessitates trimming or removal of these compromised sequences.

**12. Some quality measures are expected to give warnings for our data. Can you explain why?**

This due to the poor base quality scores that occur while sequencing and also because of the presence of adapter sequences.

**13. What is the meaning of a backslash at the end of a line in the provided Trimmomatic code?**
```bash
for file in *_1.fastq; do sample=$(basename $file _1.fastq); java -jar
../../bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -trimlog
../2_trimming/$sample.log -baseout ../2_trimming/$sample.fastq ${sample}_1.fastq
${sample}_2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:8:15
MINLEN:140 AVGQUAL:2 > ../2_trimming/${sample}_summary.log 2>&1; done
```

It is just to avoid writing all the code in one line. It helps to continue the code in separate lines.

**14. What percentage of the reads was discarded?**

```bash
/amplicon_seq/data/Colonies$ total_reads=$(for file in *.fastq; do cat "$file" | grep "^@M00" | wc -l; done | awk '{sum+=$1} END {print sum}')

echo $total_reads 
7850134

/amplicon_seq/results/02_trimmed_fastq$ discarded_reads=$(for file in *U.fastq; do cat $file | grep "^@M00" | wc -l; done | awk '{sum+=$1} END {print sum}')

echo $discarded_reads
34392

/amplicon_seq/results/02_trimmed_fastq$ echo "($discarded_reads/$total_reads)*100" | bc -l
.43810717116421197300  # 0.44 % reads were discarded after trimming
```

**15. Run the FastQC command again for the trimmed output. Choose three different categories in the fastqc report and compare the graphs the program generated. Discuss how they differ before and after the trimming.**

They do not differ much. There were very small differences in the end base pair
quality scores for the forward reads and not much difference after trimming for the
reverse reads. This does not matter because we will be using paired end reads, in
either case we will end up with good quality reads.

**16. What does pandaseq do?**

It merges two reads, forward read and the reverse read.

**17. What is the average length of the merged sequences?**

```bash
for file in *_aligned.fastq; do
    cat $file
done | awk 'NR%4==2 {total+=length($0);count++} END {print "The average length is " total/count}'
The avg length is 253.015
```

**18. Describe the rationale of the Vsearch recommended workflow.**

The goal of this workflow is to preprocess sequence reads before performing chimera
checking. The preprocessing includes:
1. Preliminary dereplication of reads within each sample (to reduce redundancy
and label reads per sample).
2. Pooling of reads from all samples into a single file.
3. Dereplicating reads across all samples (removing duplicates).
4. Preclustering reads (grouping similar sequences).

**19. What is the purpose of the initial dereplication step? In which situation is it most relevant to do this?**

Purpose:
- To reduce computational complexity by collapsing identical reads.
- To speed up later steps by working with unique sequences rather than
redundant ones.
- To retain sample identity by relabeling sequences.

Relevant Situation:
- When dealing with large datasets where sequencing depth is high.
- ​When exact sequence matches (not similar sequences) are needed for downstream analysis.

**20. Investigate the resulting file using less. What does –relabel do? What happened to the headers?**

```-relabel``` - prefixes sequence headers with basename of the file to keep track of which sample they belong to.

**21. What percentage of the reads did you discard during dereplication (not how many your kept)? Make sure you understand how the reads are represented before and after dereplication.**

```bash
/amplicon_seq/results/05_precluster$ unique_seq=$(for file in *.log.log; do cat $file | grep -oP "(\d+) unique"; done |awk '{sum+=$1} END {print sum}')
amplicon_seq/results/05_precluster$ echo $unique_seq 
580341

/amplicon_seq/results/05_precluster$ total_derep_reads=$(for file in *.derep.fasta; do cat $file | grep -oP "size=\d+" | cut -f2 -d "="; done | awk '{sum+=$1} END {print sum}')
amplicon_seq/results/05_precluster$ echo $total_derep_reads 
3888919

echo "scale=2; (($total_derep_reads - $unique_seq) / $total_derep_reads) * 100" | bc -l
85.00 # 85% of the reads were discarded
```

**22. Verify that the resulting file contains the same number of reads and nucleotides as the original files added together. How many reads are in the resulting file?**

```bash
/amplicon_seq/results/05_precluster$ unique_seq=$(for file in *.log.log; do cat $file | grep -oP "(\d+) unique"; done |awk '{sum+=$1} END {print sum}')

amplicon_seq/results/05_precluster$ echo $unique_seq 
580341

/amplicon_seq/results/05_precluster$ cat all_combined.fasta | grep "^>" | wc -l
580341 # Yes the have the same number of reads


for file in *.derep.fasta; do
    cat $file | grep -v "^>" | tr -d "\n" | wc -c;
done | awk '{sum+=$1} END {print sum}'
147103366

/amplicon_seq/results/05_precluster$ cat all_combined.fasta | grep -v "^>" | tr -d '\n' | wc -c
147103366 # Yes we have the same number of nucleotides
```

**23. What is the purpose of this dereplication step?**

This step removes duplicate sequences across all samples while keeping track of the
original abundance. By doing so:
- Redundant sequences are collapsed into unique sequences.
- Abundance data is retained (via --sizein and --sizeout).●​ A mapping file (.uc) is created, allowing us to later assign unique sequences
back to their original samples.
- Computational efficiency is improved, as fewer sequences need to be processed
in downstream steps.

**24. How many reads does the most frequent sequence from the `vsearch` output represent?**

```bash
/amplicon_seq/results/05_precluster$ cat all_derep.fasta | grep -P "size=$(cat all_derep.fasta | grep -oP "size=(\d+)"| cut -d "=" -f2 | sort -nr | uniq | head -1)"
>CT_A3;size=181549 # The most frequent output is CT_A3 and has a size of 181549
```

**25. Calculate the number of reads from each sample. This can be achieved by looping over the sample names and extracting corresponding size counts for each. Make sure they correspond to the previous counts. Write down the number of reads per sample**

```bash
# We will compare if the read matchs with the all_combined.fasta
for file in *.log.log; do 
    read_name=$(basename $file .log.log); 
    cat all_combined.fasta | grep "^>$read_name" | cut -f2 -d "=" | awk -v name="$read_name" '{sum+=$1} END {print name, sum}';
done

CT_A 914006
CT_B 703861
CT_C 366721
CT_D 285553
MA_A 1105997
MA_B 470020
MA_C 42761

for file in *.log.log; do read_name=$(basename $file .log.log);  cat all_combined.fasta | grep "^>$read_name" | cut -f2 -d "=" | awk -v name="$read_name" '{sum+=$1} END {print name, sum}'; done | cut -f2 -d " " | awk '{sum+=$1} END {print sum}'

3888919 # The added results of all the individual colonies

/amplicon_seq/results/05_precluster$ cat all_derep.fasta | grep "^>"| cut -f2 -d ";" | cut -f2 -d "=" | awk '{sum+=$1} END {print sum}'
3888919 # They correspond to the same number sequences

# In the all_derep.fasta, there reads which are same so they have been added to the most similar one in the group
/amplicon_seq/results/05_precluster$ for file in *.log.log; do read_name=$(basename $file .log.log);  cat all_derep.fasta | grep "^>$read_name" | cut -f2 -d "=" | awk -v name="$read_name" '{sum+=$1} END {print name, sum}'; done | cut -f2 -d " " | awk '{sum+=$1} END {print sum}'
3888919 # Still the same numer of total reads but per colonies the number of sequencs changes
```

**26. What is the purpose of the preclustering step?**

The preclustering step groups highly similar sequences together (≥97% identity).
This serves several purposes:
- Reduces computational load: It minimizes redundancy in the dataset, making
downstream analyses like chimera checking and OTU picking much faster.
- Maintains biological signal: While some sensitivity is lost, sequences grouped
at 97% similarity generally still represent biologically meaningful clusters.
- Creates representative sequences: The centroids represent clusters, reducing
the total number of sequences that need to be checked later.

**27. How many FASTA entries do we have after preclustering?**

```bash
cat all_preclust.fasta | grep "^>" | wc -l
48148
```

**28. What is the number of centroids in the file?**

```bash
# S - Mean the centroid, alternatively H is is for hit
amplicon_seq/results/05_precluster$ cat all_preclust.uc | awk '$1=="S"' | wc -l
48148
```

**29. What would have happened if we had run first the reference step and then the de novo?**

If we had performed reference-based chimera checking before the de novo check, we
might have inadvertently removed genuine sequences that are absent from the
reference database. This could decrease the sensitivity of the subsequent de novo
chimera detection, potentially leading to misclassifications and a loss of truly novel
sequences. Essentially, the reference check could mask or eliminate real sequences
simply because they weren't in the database, making it harder for the de novo method
to identify chimeric sequences involving those novel reads.

**30. Make sure all three files were created. Describe the content of each file.**

After running the de novo chimera checking, vsearch produces three output files:
1. all.denovo.nonchimeras.fasta
    - Contains all sequences classified as non-chimeric (valid sequences).
    - These sequences passed the chimera check and will be used for further analysis.
2. all.denovo.chimeras.fasta
    - Contains all sequences classified as chimeric.
    - These are suspected to be artifacts formed from two or more parent sequences during PCR amplification.
    - These sequences will be removed from the analysis.
3. all.denovo.uchime
    - A log file with detailed chimera classification results.
    - It contains:
        - Query sequence names.
        - Chimera status (chimeric or non-chimeric).
        - Supporting scores for chimera detection.

**31. What does de novo actually mean, and why do we call this kind of chimera checking de novo?**

De novo means "from scratch" or "without external reference." De novo chimera
detection identifies chimeric sequences without using a reference database.
Instead, it analyzes sequence abundance and similarity to detect chimeras.
Chimeric sequences tend to be low-abundance and share similarity with multiple
more abundant sequences

**32. What percentage of the reads was classified as chimeric?**

```bash
amplicon_seq/results/06_chimera$ cat all.denovo.log 
WARNING: The uchime3_denovo command does not support multithreading.
Only 1 thread used.
vsearch v2.30.4_linux_x86_64, 2003.8GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file /home/inf-21-2024/projects/amplicon_seq/results/05_precluster/all_preclust.fasta 100%
12278768 nt in 48148 seqs, min 240, max 496, avg 255
Masking 100%
Sorting by abundance 100%
Counting k-mers 100%
Detecting chimeras 100%
Found 15924 (33.1%) chimeras, 32224 (66.9%) non-chimeras,
and 0 (0.0%) borderline sequences in 48148 unique sequences.
Taking abundance information into account, this corresponds to
177746 (4.6%) chimeras, 3711173 (95.4%) non-chimeras,
and 0 (0.0%) borderline sequences in 3888919 total sequences.
```
Example to illustrate:

\>RealSeq1;size=10000 (non-chimeric, very abundant)
\>RealSeq2;size=8000 (non-chimeric, very abundant)
\>Chimera1;size=5 (chimeric, low abundance)
\>Chimera2;size=3 (chimeric, low abundance)
\>Chimera3;size=2 (chimeric, low abundance)

- By unique sequences: 3 out of 5 are chimeric (60%)
- By total reads: 10 out of 18,010 reads are chimeric (0.06%)

**4.6%** of the reads are chimeric.

**33. Inspect the database. How many sequences does it contain? Approximately how long are the sequences?**

```bash
amplicon_seq/data/chimera_database$ cat gold.fasta | grep "^>" | wc -l
5181

/amplicon_seq/data/chimera_database$ cat gold.fasta | grep -v "^>" | awk '{total+=length($0);count++} END {print "The average length is " total/count}'
The average length is 1435.68
```

**34. What does each column represent?**

| Column Name | Description |
| -------- | -------- |
| 1 | Record type (C for chimera, N for non-chimera)  |
| 2 |  Query sequence ID |
| 3 |  Score (higher score = more likely a chimera) |
| 4 |  Number of parent candidates |
| 5 |  Identity with closest parent |
| 6 |  Identity with second closest parent |
| 7 |  Chimera model score |
| 8 |  Delta score |
| 9 |  Chimera decision threshold |
| 10 |   Chimeric flag (Y = chimeric, N = non-chimeric)|
| 11  |   Number of breakpoints detected|
| 12  |   Identity of closest parent sequence|
| 13  |   Identity of second closest parent sequence|
| 14  |   Closest parent database sequence ID|
| 15  |   Second closest parent database sequence ID|
| 16  |   Length of the query sequence|
| 17  |   Query sequence abundance|
| 18  |   Chimera decision (Y or N)|
| 19 |   Extra information|

**35. How many reads were discarded (in percentage) in each of the chimera checking steps?**

```bash
/amplicon_seq/results/06_chimera$ cat all.denovo.log 
WARNING: The uchime3_denovo command does not support multithreading.
Only 1 thread used.
vsearch v2.30.4_linux_x86_64, 2003.8GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file /home/inf-21-2024/projects/amplicon_seq/results/05_precluster/all_preclust.fasta 100%
12278768 nt in 48148 seqs, min 240, max 496, avg 255
Masking 100%
Sorting by abundance 100%
Counting k-mers 100%
Detecting chimeras 100%
Found 15924 (33.1%) chimeras, 32224 (66.9%) non-chimeras,
and 0 (0.0%) borderline sequences in 48148 unique sequences.
Taking abundance information into account, this corresponds to
177746 (4.6%) chimeras, 3711173 (95.4%) non-chimeras,
and 0 (0.0%) borderline sequences in 3888919 total sequences.


/amplicon_seq/results/06_chimera$ cat all.ref.log 
vsearch v2.30.4_linux_x86_64, 2003.8GB RAM, 256 cores
https://github.com/torognes/vsearch

Reading file /home/inf-21-2024/projects/amplicon_seq/data/chimera_database/gold.fasta 100%
7438266 nt in 5181 seqs, min 1205, max 1585, avg 1436
Masking 100%
Counting k-mers 100%
Creating k-mer index 100%
Detecting chimeras 100%
Found 3842 (11.9%) chimeras, 28376 (88.1%) non-chimeras,
and 6 (0.0%) borderline sequences in 32224 unique sequences.
Taking abundance information into account, this corresponds to
14877 (0.4%) chimeras, 3696280 (99.6%) non-chimeras,
and 16 (0.0%) borderline sequences in 3711173 total sequences.
```

**39. How many non-chimeric preclusters, dereplicated reads and raw reads do we have ?**

Number of non-chimeric preclusters
```bash
/amplicon_seq/results/06_chimera$ cat all.ref.nonchimeras.fasta | grep "^>" | wc -l
28376

/amplicon_seq/results/06_chimera$ cat all.ref.nonchimeras.fasta | grep -P "size=" |  cut -f2 -d "=" | awk '{sum+=$1} END {print sum}'
3696280 # This is the count of the raw reads

# This gives you the unique sequences that survived both de novo and reference-based chimera checking.
```
Number of de-replicated reads
```bash
/amplicon_seq/results/07_clustering$ cat all.nonchimeras.derep.fasta | grep "^>" | wc -l
386235
# This counts the total reads represented by the non-chimeric dereplicated sequences.

/amplicon_seq/results/07_clustering$ cat all.nonchimeras.derep.fasta | grep -P "size=" |  cut -f2 -d "=" | awk '{sum+=$1} END {print sum}'
3696280 # This is the count of the raw reads
```
Number of raw reads
```bash
/amplicon_seq/results/07_clustering$ cat all.nonchimeras.fasta | grep -P "^>" | wc -l
471806
# This counts all the original reads that belong to non-chimeric sequences.

amplicon_seq/results/07_clustering$ cat all.nonchimeras.fasta | grep -P "size=" |  cut -f2 -d "=" | awk '{sum+=$1} END {print sum}'
3696280
```
The read count (3,696,280) remains constant during 
backtracking because we are retrieving the original entries with 
their original ;size= abundance tags intact. We are not 
recalculating or duplicating reads - we are simply expanding from 
a collapsed representation (preclusters) back to more granular 
representations (dereplicated, then sample-level), while preserving 
the abundance information tracked by the size tags.

**40. For each step, compare the percentage of the total number of reads, dereplicated reads and preclusters that are classified as chimeric.**

```bash



```

**41. Study Linux sed command. How does it change the sample headers?**

sed "s/[0-9]\+;/;/" all.nonchimeras.fasta > all.nonchimeras.renamed.fasta
It just replaces the header of every sequence that starts with a digit one or more times followed by a ; with just a ;.

**42. Calculate the total number of OTUs and how many sequences they represent together. Does this correspond to the number of input sequences?**

```bash
/amplicon_seq/results/07_clustering$ cat outs.fasta | grep "^>" | wc -l
55406 # There are 55406 OTUs

/amplicon_seq/results/07_clustering$ cat outs.fasta | grep "size=" | cut -f2 -d "=" | awk '{sum+=$1} END {print sum}'
3696280 # This is the total counts of the OTUs (how many sequences they represent together)
```
Yes. It corresponds to the number of input sequences

**43. How many OTUs were identified in each of the seven colonies? (Hint: investigate the otus.tsv file)**



**44. What was the average number of sequences per OTU in each of the seven colonies?**


**45. Explain what the RDP database is.**

The RDP is a curated database of aligned and annotated rRNA gene sequences
(mainly 16S rRNA) used for microbial ecology and taxonomy. It provides tools for
sequence classification, alignment, and phylogenetic analysis.

**46. How many sequences are in it?**

The number of sequences varies by version. The RDP database as of 27th Nov 2013
contains 2,809,406 aligned and annotated bacterial and archaeal small subunit (SSU)
rRNA gene sequences and 62,860 fungal large subunit (LSU) rRNA gene sequences.


**47. Explain how can the RDP database increase or decrease in size in future releases.**

The database size changes due to addition of new sequences, improved
taxonomy/alignment, removal of low-quality sequences, user submissions, and
taxonomic classification changes.

**48. How many cluster sequences could not be classified as bacterial and be determined on the phylum level (on the 0.8 confidence level)?**

```bash
/amplicon_seq/results/07_clustering$ cat all.fixedRank | cut -f3 | sort -n | uniq -c
     57 Archaea
  55346 Bacteria
      3 Eukaryota
# 60 were not classified as bacteria

/amplicon_seq/results/07_clustering$ cat all.fixedRank | awk -F "\t" '{sum+=$8>=0.8} END {print sum}'
37284
(base) inf-21-2024@bioinf-serv2:~/projects/amplicon_seq/results/07_clustering$ cat all.fixedRank | awk -F "\t" '$8>=0.8' | wc -l
37284
```

**49. How many cluster sequences were determined on the genus level (using 0.8)?**

```bash
/amplicon_seq/results/07_clustering$ cat all.fixedRank | awk -F "\t" '$20>=0.8' | wc -l
9354
```

**50. How would it change if you had used thresholds of 0.7 and 0.9?**

If the threshold is 0.7, the number of cluster sequences would increase and if the threshold is 0.9, the number of cluster sequences would decreases with respect to 0.8 threshold.
