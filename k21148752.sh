#!/bin/bash #

## Here we are making the new directories to run the pipeline. The directories are called ngs_course, dnaseq, data, meta, results and, log. The log directory can be useful to keep track of the commands that were used.
mkdir ngs_course
mkdir ngs_course/dnaseq
cd ngs_course/dnaseq
mkdir data meta results log
ls -lF
cd ~/ngs_course/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq

## Here the raw fastq data used for the pipeline is being downloaded.
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

## The fastq files are then copied into the untrimmed_fastq folder and the bed file into the data directory.
mv *fastq.qz ~/ngs_course/dnaseq/data/untrimmed_fastq
mv annotation.bed ~/ngs_course/dnaseq/data

## Here the reference files which we will map against in the alignment step are downloaded and moved into the data folder.
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ~/ngs_course/dnaseq/data/

## Here the tools requiured for to run the pipeline are downloaded.
cd ~/
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./ Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib

## Here we are reloacting to the untrimmed fastq folder and verifying it's content.
cd ~/ngs_course/dnaseq/data/untrimmed_fastq
ls -lart

## This command, the .qz files are decompressed and then combined. After this is done, fastqc step is done for the quality assessment of the raw reads.
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq
fastqc NGS0001.R1.fastq NGS0001.R2.fastq

## A new directory is made called fastqc untrimmed reads, it is here that the fastqc results will be stored
mkdir ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/

cd ~/ngs_course/dnaseq/data/untrimmed_fastq

## The trimmomatic step is undertaken to trim the adapters from the raw data
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq \
  -baseout /home/ubuntu/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R \
 ILLUMINACLIP:/home/ubuntu//miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

## Quality assessment is undertaken on the trimmed reads
fastqc ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P
fastqc ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_2P

## A new directory by the name of reference is created
mkdir -p ~/ngs_course/dnaseq/data/reference
mv ~/ngs_course/dnaseq/data/hg19.fa.gz ~/ngs_course/dnaseq/data/reference/

## Here bwa indexes sequences in the FASTA format
bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz
ls ~/ngs_course/dnaseq/data/reference

## The untrimmed_fastq folder is no longer needed and is therefore deleted to create more space on the virtual machine.
rm -r ~/ngs_course/dnaseq/data/untrimmed_fastq

mkdir ~/ngs_course/dnaseq/data/aligned_data

## BWA-MEM is an alignment algorithm for aligning sequence reads against a reference genome

bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_2P > ~/ngs_course/dnaseq/data/aligned_data/NGS0001.sam

cd ~/ngs_course/dnaseq/data/aligned_data

## Here the sam file previously downloaded is being converted into a .bam (binary) version
samtools view -h -b NGS0001.sam > NGS0001.bam

## The .bam file is then sorted
samtools sort NGS0001.bam > NGS0001_sorted.bam

## The sorted bam file is then indexed
samtools index NGS0001_sorted.bam


cd ~/ngs_course/dnaseq/data/aligned_data

## Flagstat and idxstats are alignment statistics tools which produce certain types of statistical data.
samtools flagstat NGS0001_sorted.bam
samtools idxstats NGS0001_sorted.bam

## The NGS0001.sam file is a large file that is no longer needed and can therefore be deleted to free up more space on the virtual machine.
cd ~/ngs_course/dnaseq/data/aligned_data
rm NGS0001.sam

## Picard is then used to mark duplicate reads in the sorted.bam file
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

## here the sorted_marked.bam file is indexed
samtools index NGS0001_sorted_marked.bam

## The reads are then filtered using a set criteria in which the Minimum MAPQ quality score of 20 is set and a bitwise flag of 1796. 
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

## here the sorted_filtered.bam file is indexed
samtools index NGS0001_sorted_filtered.bam

## The reference genome is decompressed
zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa

## The reference genome is being indexed using samtools faidx
samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa

## Freebayes is used to call variants
freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/NGS0001.vcf

## bgzip is used to compress the vcf file
bgzip ~/ngs_course/dnaseq/results/NGS0001.vcf

## tabix is used to index file
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001.vcf.gz

## The vcf file is then filtered using their quality scores as an indicator
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf

cd ~/ngs_course/dnaseq/data/aligned_data

## bedtools intersect has the role of screening the overlaps between the vcf file and annotation bed
bedtools intersect -header -wa -a ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf -b ../annotation.bed > ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

## bgzip is used to compress the filtered_R vcf file
bgzip ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

## tabix is used to index the filtered_R vcf file
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf.gz

## The trimmed_fastq folder is no longer required and can therefore be deleted to create more space on the virtual machine
rm -r ~/ngs_course/dnaseq/data/trimmed_fastq

## Annovar is downloaded using the annovar website and then transferred over to the home folder on the virtual machine using FileZilla.
cd ~/

## The annovar file is then unzipped
tar -zxvf annovar.latest.tar.gz
cd annovar
## Numerous different databases of annovar are then downloaded
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/

## The VCF file is then converted into a annovar input file, in the process it ensures that any VCF specific information is not lost
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_R.avinput

## table annovar uses the VCF file to generate a tab-delimited file with numerous sets of data, in this particular it is the CSV file that is produced
./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_R.avinput humandb/ -buildver hg19 -out ~/ngs_course/dnaseq/results/NGS0001 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

