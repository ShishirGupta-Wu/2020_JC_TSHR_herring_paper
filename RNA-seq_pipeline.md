# Pipeline from QC of raw reads to assignment of read counts to genes

## 1. QC of raw reads with Trimmomatic

The new script `02-trim_raw_reads-SE.sh` looks like this:
```bash
#!/bin/bash
#SBATCH -A XXX
#SBATCH -M snowy
#SBATCH -p core -n 6
#SBATCH -t 1-00:00:00
#SBATCH -J arr_trim_RNA_spawn
#SBATCH -e trim_RNA_spawn_%J_%A_%a.err   
#SBATCH -o trim_RNA_spawn_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-16

# Load libraries.
module load bioinfo-tools
module load trimmomatic/0.36
module load FastQC/0.11.8

echo This is array job: $SLURM_ARRAY_TASK_ID

R1_FILE='./R1_files.txt'
ID_FILE='./sample_IDs.txt'
OUT_DIR='/PATH/Transcriptome_Data/Spawning/clean-reads'
OUT_SUBDIR='./sample_subDir.txt'
ADAPTERS_FILE='/PATH/Transcriptome_Data/Spawning/clean-reads/adapters.fa'

R1_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $R1_FILE)
ID_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_FILE)
SUBDIR_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $OUT_SUBDIR)

echo -e This is the R1 target: ${R1_target}
echo -e This is the ID target: ${ID_target}
echo -e This is the SUBDIR target: ${OUT_DIR}/${SUBDIR_target}
echo -e This is the R1 output file name: ${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R1.fq.gz

# Create a sub-folder for each sample
mkdir ${OUT_DIR}/${SUBDIR_target}

# Run trimmomatic on raw RNA sequence reads - Single-End
java -Xmx38g -jar /sw/apps/bioinfo/trimmomatic/0.36/rackham/trimmomatic-0.36.jar SE -threads 6 -phred33 \
${R1_target} \
${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R1.fq.gz \
ILLUMINACLIP:${ADAPTERS_FILE}:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:15 TRAILING:15 MINLEN:36

# Run fastQC on the trimmed reads file
fastqc ${OUT_DIR}/${SUBDIR_target}/${ID_target}_trimmed_R1.fq.gz -o ${OUT_DIR}/FastQC_results

```
Next, I obtained a single FastQC report for all samples using MultiQC:
```
module load bioinfo-tools
module load MultiQC/1.7

cd /PATH/Transcriptome_Data/Spawning/clean-reads/FastQC_results

multiqc .
```

## 2. Mapping of reads to the genome
### 2.1. Generate genome index file

Code used for the standard herring genome (chr + unplaced scaffolds), called `03-1-generate-genomeIndex-gsnap_A.sh`:
```bash
#!/bin/bash
#SBATCH -A XXX
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J genomeIndex_gsnapA
#SBATCH -e genomeIndex_gsnapA_%J_%A_%a.err
#SBATCH -o genomeIndex_gsnapA_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL

# Set environment variables.
# Files.
WORK_DIR='/PATH/RNA-seq/genome/gsnap'
GENOME_NAME='Ch_v2.0.2.fasta'
GENOME_PATH='/PATH/RNA-seq/genome'

# Load required programs.
module load gmap-gsnap/2018-07-04 samtools

# Create a genome index. Default parameters, except -s none (turn sorting off), -D (set output dir)
gmap_build -d $GENOME_NAME.gsnap -D $WORK_DIR $GENOME_PATH/$GENOME_NAME -s none

```

### 2.3 Create an IIT file with SNP information for conducting a SNP-tolerant alignment

For the standard Herring genome, used `03-2-generate-SNPindex-gsnap_A.sh`:
```bash
#!/bin/bash
#SBATCH -A XXX
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10:00:00
#SBATCH -J SNPindex_gsnapA
#SBATCH -e SNPindex_gsnapA_%J_%A_%a.err
#SBATCH -o SNPindex_gsnapA_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL

# Load required programs.
module load bioinfo-tools
module load gmap-gsnap/2018-07-04
module load samtools/1.10

# Set environment variables.
# Files.
GSNAP_GENOME_PATH='/PATH/RNA-seq/genome/gsnap/Ch_v2.0.2.fasta.gsnap'
GSNAP_GENOME_NAME='Ch_v2.0.2.fasta.gsnap'
VCF_FILE='/PATH/SNP_populations/67.SNP.vcf'
VCF_NAME='pool_67.hf.SNPs'

# Go to working directory.
cd $GSNAP_GENOME_PATH/*.maps/

# Create a txt file with SNP info from a VCF file.
vcf_iit $VCF_FILE > $VCF_NAME.txt

# Convert the txt file into a IIT file.
cat $VCF_NAME.txt | iit_store -o $VCF_NAME

# Create a reference space index and compressed genome.
snpindex -d $GSNAP_GENOME_NAME -D $GSNAP_GENOME_PATH -V $GSNAP_GENOME_PATH/*.maps -v $VCF_NAME.SNPsdb $VCF_NAME.iit

```


### 2.4 Map RNA-seq reads with GSNAP

- The script used for the standard Herring genome is called `03-3-readMapping-gsnap-snpT-SE_A.sh`:
```bash
#!/bin/bash
#SBATCH -A XXX
#SBATCH -p node
#SBATCH -t 02-00:00:00
#SBATCH -J readMap_gsnapA
#SBATCH -e readMap_gsnapA_%J_%A_%a.err
#SBATCH -o readMap_gsnapA_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -a 1-16

# Load required programs.
module load bioinfo-tools
module load gmap-gsnap/2018-07-04
module load samtools/1.10

# Set environment variables.
# Files.
WORK_DIR='/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap'
GSNAP_GENOME_PATH='/PATH/RNA-seq/genome/gsnap/Ch_v2.0.2.fasta.gsnap'
GSNAP_GENOME_NAME='Ch_v2.0.2.fasta.gsnap'
VCF_NAME='pool_67.hf.SNPs'
FASTQ_LIST='/PATH/RNA-seq/TSHR/bams/fastq_files_BSH-BR.list'
ID_LIST='/PATH/RNA-seq/TSHR/bams/sample_IDs_BSH-BR.list'

# For a given job in the array, set the correspondent sample files.
FASTQ_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FASTQ_LIST)
ID_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_LIST)

# Print current file info to stdout for future reference.
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FASTQ_target: ${FASTQ_target}
echo -e This is the ID_target: ${ID_target}

# Create directory for current sample and temporal files.
if [ -d $WORK_DIR/${ID_target} ]; then echo "tmp/ exists in WORK_DIR"; else mkdir $WORK_DIR/${ID_target}; fi
# Create directory for temporal files.
if [ -d $WORK_DIR/${ID_target}/tmp ]; then echo "tmp/ exists in WORK_DIR"; else mkdir $WORK_DIR/${ID_target}/tmp; fi

# Directory where the BAM files will be stored.
cd $WORK_DIR/${ID_target}

# Map RNA-seq reads with GSNAP
gsnap \
--gunzip \
-m 0.1 \
-D $GSNAP_GENOME_PATH \
-d $GSNAP_GENOME_NAME \
-V $GSNAP_GENOME_PATH/*.maps \
-v $VCF_NAME.SNPsdb \
-B 5 -N 1 -n 30 -E 4 \
--nthreads=16 \
--gmap-mode=pairsearch,ends,improve \
-A sam -J 33 -O --quiet-if-excessive \
--read-group-id=${ID_target} \
--read-group-name=${ID_target} \
--read-group-library=50SE \
--read-group-platform=Illumina ${FASTQ_target} | samtools view -@ 16 -F 4 -q 20 -Sb - | samtools sort -@ 16 -T $WORK_DIR/${ID_target}/tmp - > ${ID_target}.gsnap.snpT.sort.bam

# Create an index for the bam file
samtools index ${ID_target}.gsnap.snpT.sort.bam
samtools stats ${ID_target}.gsnap.snpT.sort.bam > ${ID_target}.stat

```

## 3. Assign read counts to genes

Script `04-reads2genes-fC-SE_withoutPrimary.sh`:
```bash
#!/bin/bash
#SBATCH -A XXX
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 5:00:00
#SBATCH -J reads2genes_A
#SBATCH -e reads2genes_A_%J_%A_%a.err
#SBATCH -o reads2genes_A_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL

# Load required programs.
module load bioinfo-tools
module load subread/2.0.0

# Set environment variables.
# Files.
WORK_DIR='/PATH/RNA-seq/TSHR/analysis/04-featureCounts'
OUTPUT_NAME='Her-16samples-RNAseq-BSH-BR'
GTF_FILE='/PATH/RNA-seq/genome/Clupea_harengus.Ch_v2.0.2.96.renamedChr.gtf'

# Go to the working directory
cd $WORK_DIR

# Assign reads to genes using featureCounts.
featureCounts \
-s 2 -t exon -g gene_id \
-T 2 -Q 20 \
-o $OUTPUT_NAME.gsnap.s2.q20.gene_id.featureCounts.txt \
-a $GTF_FILE \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T10F-2-BR/Sample_Ber16-T10F-2-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T10F-BSH/Sample_Ber16-T10F-BSH.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T13F-3-BR/Sample_Ber16-T13F-3-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T13F-BSH/Sample_Ber16-T13F-BSH.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T17F-BR/Sample_Ber16-T17F-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T17F-BSH/Sample_Ber16-T17F-BSH.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T18F-BR/Sample_Ber16-T18F-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_Ber16-T18F-BSH/Sample_Ber16-T18F-BSH.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-14F-BR/Sample_HAS-14F-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-14F-BSH/Sample_HAS-14F-BSH.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-5F-BR/Sample_HAS-5F-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-5F-BSH/Sample_HAS-5F-BSH.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-8F-BR/Sample_HAS-8F-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-8F-BSH/Sample_HAS-8F-BSH.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-9F-BR/Sample_HAS-9F-BR.gsnap.snpT.sort.bam \
/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_HAS-9F-BSH/Sample_HAS-9F-BSH.gsnap.snpT.sort.bam

#/PATH/RNA-seq/TSHR/bams/Ch_v2.0.2.fasta.gsnap/Sample_*/*.bam

```
