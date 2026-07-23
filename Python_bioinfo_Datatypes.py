##strings

##manipulation of strings

##1. replace(old,new[,count])
dna_seq = "GGTTGCGCTAGCTCTGA"
rna_seq = dna_seq.replace('T', 'U')
print(rna_seq)


##2. count(sub[,start[,end]])
## can be used to calculate CG content
c = dna_seq.count("C")
g = dna_seq.count("G")
print((c+g)/len(dna_seq)*100)

##3. find(sub[,start[,end]])
##returns position of the substring between start and end positions else returns -1
print(rna_seq.find("CGA"))
print(dna_seq.find("AAG"))

##4. index(sub[,start[,end]])
##similar to find function but raises a ValueError exception if absent. 
print(rna_seq.index("CGA"))
print(dna_seq.index("AAG"))

##5. split
##separates words of a string and returns a list default is white space 
print('dna_seq'.split())
print('rna_seq'.split(","))

##6. join
##reverse of split function 


##B. LISTS
##collection of ordered objects enclosed in square brackets []

##1. 
sample_list = [1, "two", 4, "seven"]
print("sample_list")
##2. 
next_list = ["apple", "banana"]
print("next_list")
##3. 
nested_list = [next_list, sample_list]
print("nested_list")

##4. accessing list elements
nested_list = nested_list[3]
sample_list = sample_list[-1]

##5. converting non-list object to list

aseq = "GGAATGTGCGAGTGTGC"
list(aseq)

##list comprehension 

##addition
##1. append(element)
##2. insert(position, element)
##3. extend(list) - similar to + 

##remove
##1. pop(index)
##2. remove(element)
##3. del

##copying list
##1. copy() from the copy module
import copy 

##2. 