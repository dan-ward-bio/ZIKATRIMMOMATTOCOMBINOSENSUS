#trim primers off  reads using trimmomatic

trimmomatic.jar PE -phred33 -trimlog trimlog.txt POOL1_R1.fastq.gz POOL1_R2.fastq.gz POOL1_TRIMMED_R1_PAIRED.fq.gz POOL1_TRIMMED_R1_UNPAIRED.fq.gz POOL1_TRIMMED_R2_PAIRED.fq.gz POOL1_TRIMMED_R2_UNPAIRED.fq.gz ILLUMINACLIP:./primers1.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic.jar PE -phred33 -trimlog trimlog.txt POOL2_R1.fastq.gz POOL2_R2.fastq.gz POOL2_TRIMMED_R1_PAIRED.fq.gz POOL2_TRIMMED_R1_UNPAIRED.fq.gz POOL2_TRIMMED_R2_PAIRED.fq.gz POOL2_TRIMMED_R2_UNPAIRED.fq.gz ILLUMINACLIP:./primers2.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzip -d POOL1_TRIMMED_R1_PAIRED.fq.gz POOL1_TRIMMED_R2_PAIRED.fq.gz POOL2_TRIMMED_R1_PAIRED.fq.gz POOL2_TRIMMED_R2_PAIRED.fq.gz

cat POOL1_TRIMMED_R1_PAIRED.fq POOL2_TRIMMED_R1_PAIRED.fq > COMBINED_TRIMMED_R1.fq

cat POOL1_TRIMMED_R2_PAIRED.fq POOL2_TRIMMED_R2_PAIRED.fq > COMBINED_TRIMMED_R2.fq

gzip COMBINED_TRIMMED_R1.fq

gzip COMBINED_TRIMMED_R2.fq

#Index ref ensuring there are no gaps in the sequence and that the sequence is appropriate, taking in to accoun things like the presence or lack of the 5' and 3' UTR

bwa index ref.fasta

#Map using bwa mem

bwa mem ref.fasta COMBINED_TRIMMED_R1.fq.gz COMBINED_TRIMMED_R2.fq.gz > MAPPED.sam

#convert SAM to BAM

samtools view -q 15 -b -S MAPPED.sam > MAPPED.bam

#sort the BAM

samtools sort MAPPED.bam -o SORTED.bam

#index the sorted bam vor visulisation in Tablet

samtools index SORTED.bam

#now you can open the sorted bam and the indexed reference in tablet

#BEDTOOLS to find where there is zero coverage

samtools faidx  ref.fasta

# R-COVERAGE-OUTPUT-bedtools genomecov -dz -ibam SORTED.bam -g ref.fasta.fai > GENOME_COVERAGE_GRAPH_DATA


bedtools genomecov -bga -ibam SORTED.bam -g ref.fasta.fai > bedtools_output.bed

#Extract the positions where there is zero coverage

cat bedtools_output.bed | awk '$4 < 100' > ZERO_COVERAGE_POSITIONS.bed

#now we want to do the SNP calling with mpileup (I think the ref needs to be faidx indexed)

bcftools mpileup -Ou -f ref.fasta SORTED.bam | bcftools call -mv -Oz -o calls.vcf.gz 

#then run this command

tabix calls.vcf.gz

#index the SNP calls

bcftools index calls.vcf.gz 

#now create the consensus sequence

bcftools consensus -f ref.fasta -m ZERO_COVERAGE_POSITIONS.bed  calls.vcf.gz -o CONSENSUS_GAPS.fasta

#Tidy up

rm bedtools_output.bed calls.vcf* MAPPED* ref.fasta.* ZERO_COVERAGE_POSITIONS.bed  













