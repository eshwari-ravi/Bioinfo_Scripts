##awk,grep,sed pattern searching

##using awk for file manipulation


### using sed to replace 


###using grep to 



##parallelization using shell

##running scripts in server  in parallel 



##automation 


##checking disk storage file size 


##searching for files 


##commands with filters for bioinfo pipelines 


##shell scripting 

##1. awk, grep, sed, wc

##to fetch the count in a file
wc -l 
##to fetch total no. of sequences 
grep ">" fastafile.fa | wc -l

OR

grep -c ">" fastafile.fa

##to extract header sequences from fastq file

grep "^@" fastqfile.fq | awk '{print $1}'

##renaming files can be done using
mv "oldfilename.fa" "newfilename.fa"

## to locate a program file in the users path 
which conda

which -v
which --version 

##to find the version of a tool installed 

conda --version 






##using seqkit for fasta/fastq file manipulations
Seqkit is used for 1. BAM file processing, 2. Multiple sum or merging, 3. ordering i.e. sort, shuffle,
4. searching by grep , locate, 5. Editing - replace, rename, concat, restart, mutate, sana 6. format conversion - fa2fq, convert, tab2fx, tab2fa, fq2fa
7. basics - subseq, faidx, seq, stats, sliding, scat, translate, scat, watch
8. set operation - split, head, sample, common, rmdup, 

##

