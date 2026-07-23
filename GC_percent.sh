## calculate gc content in a fasta file using shell 

#!bin/bash

##ask user to enter a sequence 

echo "Please enter a DNA sequence"

read DNA 

##calculate total G and C in the seq

G = $(echo "$DNA" | grep -o "G" | wc -l)
C = $(echo "$DNA" | grep -o "C" | wc -l) 
##print the total length of the input seq
total = $(#DNA)

GC_percent = (echo "scale=2; (($G+$C) / $total) * 100" | bc)

##print the output
echo "total gc content in the input sequence is : $GC_percent"




