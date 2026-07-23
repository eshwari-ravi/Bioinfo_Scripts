#!/bin/bash

echo "please enter your name"

read name


echo "Hi $name, please welcome to our website"


## in python

print(f"Hi, please enter your name")
name = input()
print(f"Hello, {name}")



###function for generating gc_percent in shell 
gc_percent() {
    local seq=$1
    seq=$(echo "$seq | tr 'a-z' 'A-Z)
    local g=$(echo $seq | grep -o "G" | wc -l)
    local c=$(echo $seq | grep -o "C" | wc -l)
    local total=${#seq}
    awk -v g="$g" -v c="$c" -v t="$total" 'BEGIN {printf"%.2f", ((g+c)/t*100)}'
}
gc=$(gc_percent"AATTGGCGCGC")
echo "GC Content $gc"


#!/bin/bash

# Set input and output directories
INPUT_DIR="data"
OUTPUT_DIR="results"
REFERENCE="reference.fa"  # Reference genome

mkdir -p "$OUTPUT_DIR"

# Loop through all *_R1.fastq files
for R1 in "$INPUT_DIR"/*_R1*.fastq; do
    # Get base sample name (remove _R1.fastq)
    BASE=$(basename "$R1" | sed 's/_R1.*\.fastq//')

    # Find matching R2 file
    R2="$INPUT_DIR/${BASE}_R2.fastq"
    
    if [[ -f "$R2" ]]; then
        echo "Processing sample: $BASE"
        
        # Define output
        SAM_OUT="$OUTPUT_DIR/${BASE}.sam"
        
        # Example command: BWA MEM alignment
        bwa mem "$REFERENCE" "$R1" "$R2" > "$SAM_OUT"
        
        # Add your other steps (e.g., SAM to BAM, GATK, etc.)
        
    else
        echo "Missing R2 file for $BASE, skipping."
    fi
done



# Convert SAM to sorted BAM
BAM="$OUTPUT_DIR/${BASE}.sorted.bam"
samtools view -bS "$SAM_OUT" | samtools sort -o "$BAM"

# Index BAM
samtools index "$BAM"

# Call variants using GATK
VCF="$OUTPUT_DIR/${BASE}.vcf"
gatk HaplotypeCaller -R "$REFERENCE" -I "$BAM" -O "$VCF"

# Annotate variants using SnpEff
ANN_VCF="$OUTPUT_DIR/${BASE}.annotated.vcf"
snpEff ann your_db "$VCF" > "$ANN_VCF"

