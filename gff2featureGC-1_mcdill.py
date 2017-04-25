#!/usr/bin/env python3.6

import sys


#a function to clean up a DNA sequence
def clean_seq(input_seq):
	clean = input_seq.upper()
	clean = clean.replace('N', '')
	return clean
	
def nuc_freq(sequence, base, sig_digs=2):
	#calculate length of sequence
	length = len(sequence)
	
	#count the number of this nucleotide
	count_of_base = sequence.count(base)
	
	#calculate the base frequency
	freq_of_base = count_of_base/length 
	
	#return frequency and length 
	return(length, round(freq_of_base, sig_digs))
	
#for above: sig_digs equal to 2
#gc content for each category	

#Thoughts: cat both files (or hard code), define both with variable, read separately,
#count feature types, count Cs and Gs total, break it down by feature type, worry about formatting

#declare fsa file name
fsa_filename ='watermelon.fsa'
fsa_infile = open(fsa_filename, 'r')

#declare variable
genome = ''

line_number = 0

#read in genome file
for line in fsa_infile:
	#print(str(line_number) + ": " + line)
	
	#remove newline
	line = line.rstrip('\n')
	
	if line_number > 0:
		genome = genome + line 
		
	#increment line number
	line_number += 1
	
#print(len(genome))	
#seems right (same as above)

fsa_infile.close()

#declare gff file name
gff_filename = 'watermelon.gff'
gff_infile = open(gff_filename, 'r')

cds = ''
trna = ''
rrna = ''
intron = ''
misc = ''
repeats = ''

#read in gff file
for line in gff_infile:

	line = line.rstrip('\n')
	
	types = line.split('; type ')
	other_type = types[len(types)-1]
	#print(other_type)

	fields = line.split('\t')
	type = fields[2]
	start = int(fields[3])
	stop = int(fields[4])
	#print(type, start, "\t", stop)
	
	fragment = genome[start-1:stop]
	
	fragment = clean_seq(fragment) 
	#print(clean)
	#sys.exit()
	
	#print(fragment)
	
	if type == 'CDS':
		cds += fragment
		#cds equals big concatenated cds, etc
		
	if type == 'intron':
		intron += fragment
		
	if type == 'misc_feature':
		misc += fragment
		
	if type == 'repeat_region':
		repeats += fragment
		
	if type == 'rRNA':
		rrna += fragment
		
	if type== 'tRNA':
		trna += fragment
		
list_of_features = ['cds', 'intron', 'misc', 'repeats', 'rrna', 'trna']
feature_sequences = [cds, intron, misc, repeats, rrna, trna]
bases = ['A','G','T','C']

for i in range(len(list_of_features)):
	print(i)
	#loop over four nucleotides		
	for nucleotide in bases:
	
		#calculate the nucleotide composition for each feature
		(feature_length, feature_comp) = nuc_freq(feature_sequences[i], base=nucleotide, sig_digs=2)
		print(list_of_features[i] + "\t" + str(feature_length) + "\t" + str(feature_comp) + str(nucleotide))

##print out nucleotide frequency
##ch08
		
#and print
#print(cds.count('G'))
#print(cds.count('C'))

GC_content_cds = (cds.count('C') + cds.count('G'))
#print(GC_content_cds)

GC_content_intron = (cds.count('C') + intron.count('G'))
#print(GC_content_intron)

GC_content_misc = (misc.count('C') + misc.count('G'))
#print(GC_content_misc)

GC_content_repeats = (repeats.count('C') + repeats.count('G'))
#print(GC_content_repeats)

GC_content_rrna = (rrna.count('C') + rrna.count('G'))
#print(GC_content_rrna)

GC_content_trna = (trna.count('C') + trna.count('G'))
#print(GC_content_trna)





#print("cds \t" + str(GC_content_cds) + " \t" + str(GC_content_cds/len(genome)))
#print("intron \t" + str(GC_content_intron) + " \t" + str(GC_content_intron/len(genome)))
#print("misc_feature \t" + str(GC_content_misc) + " \t" + str(GC_content_misc/len(genome)))
#print("repeat_region \t" + str(GC_content_repeats) + " \t" + str(GC_content_repeats/len(genome)))
#print("rRNA \t" + str(GC_content_rrna) + " \t" + str(GC_content_rrna/len(genome)))
#print("trna \t" + str(GC_content_trna) + " \t" + str(GC_content_trna/len(genome)))

#"{:%.1f}%.format(GC_content_cds)/(len(genome))")

							
	

gff_infile.close()
	
		
		
		
		
	
	
