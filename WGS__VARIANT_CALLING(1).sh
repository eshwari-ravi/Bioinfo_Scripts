mkdir 1_RawReads 2_TrimmedReads 3_Alignment 4_Variant_calling
mv *.gz 1_RawReads
cd 1_RawReads
for i in `ls *_R1.fastq.gz|cut -f1 -d "_"`
do
fastp -a AGATCGGAAGAGC -w 16 -x -g -i $i"_R1.fastq.gz" -I $i"_R2.fastq.gz" -o "../2_TrimmedReads/"$i"_Trimmed_R1.fq.gz"  -O "../2_TrimmedReads/"$i"_Trimmed_R2.fq.gz"
done

bwa index Staphylococcus_aureus_ATCC_25923.fa
bwa mem -M -t 20 Staphylococcus_aureus_ATCC_25923.fa Sample2-8-12-21-A_Trimmed_R1.fq.gz Sample2-8-12-21-A_Trimmed_R2.fq.gz| samtools view -@ 20 -bS >Sample2-8-12-21-A.bam
samtools sort -@ 20 Sample2-8-12-21-A.bam -o Sample2-8-12-21-A.sorted.bam


sambamba markdup -r -t 20 LBB.sorted.bam LBB_rmdup.bam

for i in `ls *.sorted.bam| cut -f 1 -d "."`
do
freebayes -f Staphylococcus_aureus_ATCC_25923.fa -C 2 -F 0.05 -0 -q 20 -p 1 $i".sorted.bam" | vcffilter -f "QUAL > 20" >$i".vcf"
done


java -jar snpEff.jar build -gff3 -v <db name>
java -Xmx8g -jar ../RIDHI/snpEff/snpEff.jar -c ../RIDHI/snpEff/snpEff.config C.auris Sample-1-Treated.vcf> Sample-1-Treated_snpeff_annot.vcf
java -jar ../RIDHI/snpEff/SnpSift.jar extractFields Sample-1-Treated_snpeff_annot.vcf CHROM POS ID REF ALT FILTER AF AC DP MQ "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].IMPACT:" "ANN[0].GENE:" "ANN[0].GENEID">Sample-1-Treated_snpeff_annot_SIFT.vcf
