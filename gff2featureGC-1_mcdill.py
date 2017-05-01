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
	
	# key = feature type, value = concatenation of all sequences of that type
	#not useful for anything other than calculating AT/GC content
feature_sequences = {}

# key = gene name, value = another dictionary [key = exon number, value = exon sequence]
# gene_sequences[cox1][1] = 'the sequence for the first exon of cox1'
# gene_sequences[cox1][2] = 'the sequence for the second exon of cox1'
gene_sequences = {}
	
#for above: sig_digs equal to 2
#gc content for each category	


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
	
#extract and clean the sequence of this feature from the genome
	
	fragment = genome[start-1:stop]
	fragment = clean_seq(fragment) 
	#print(clean)
	#sys.exit()
	
	# determine the strand, reverse complement or not
	if(fields[6] == '-'):
		fragment = reverse_complement(fragment)
	
	#store the big concatenated thing for calculating GC content
	
	if type in feature_sequences:
		feature_sequences[type] += fragment
	else:
		feature_sequences[type] = fragment
		


if(type == 'CDS'):	
	#get the gene name
	#print(fields[8])
	attributes = fields[8].split(';')
	print(attributes[0])
	
	gene_fields = attributes[0].split(' ')
	gene_name = gene_fields[1]
	
	# get the exon#
	if('exon' in gene_fields ):
		#print("has exons: " + attributes[0])
		exon_num = gene_fields[-1]
		print(gene_name, exon_num)
	# else:
		# print("Doesn't have exons: " + attributes[0])
			
	
	#print(type, "\t", start, "\t", end)
	
		
#close the gff file	
gff_infile.close()

for feature, sequence in feature_sequences.items():
	print(feature + "\t" + str(len(sequence)))


	#print(fragment)
	
	#if type == 'CDS':
		#cds += fragment
		#cds equals big concatenated cds, etc
		
	#if type == 'intron':
		#intron += fragment
		
	#if type == 'misc_feature':
		#misc += fragment
		
	#if type == 'repeat_region':
		#repeats += fragment
		
	#if type == 'rRNA':
		#rrna += fragment
		
	#if type== 'tRNA':
		#trna += fragment
		
#list_of_features = ['cds', 'intron', 'misc', 'repeats', 'rrna', 'trna']
#feature_sequences = [cds, intron, misc, repeats, rrna, trna]
#bases = ['A','G','T','C']

#for i in range(len(list_of_features)):
	#print(i)
	#loop over four nucleotides		
	#for nucleotide in bases:
	
		#calculate the nucleotide composition for each feature
		#(feature_length, feature_comp) = nuc_freq(feature_sequences[i], base=nucleotide, sig_digs=2)
		#print(list_of_features[i] + "\t" + str(feature_length) + "\t" + str(feature_comp) + str(nucleotide))

##print out nucleotide frequency
##ch08
		
#and print
#print(cds.count('G'))
#print(cds.count('C'))

#GC_content_cds = (cds.count('C') + cds.count('G'))
#print(GC_content_cds)

#GC_content_intron = (cds.count('C') + intron.count('G'))
#print(GC_content_intron)

#GC_content_misc = (misc.count('C') + misc.count('G'))
#print(GC_content_misc)

#GC_content_repeats = (repeats.count('C') + repeats.count('G'))
#print(GC_content_repeats)

#GC_content_rrna = (rrna.count('C') + rrna.count('G'))
#print(GC_content_rrna)

#GC_content_trna = (trna.count('C') + trna.count('G'))
#print(GC_content_trna)





#print("cds \t" + str(GC_content_cds) + " \t" + str(GC_content_cds/len(genome)))
#print("intron \t" + str(GC_content_intron) + " \t" + str(GC_content_intron/len(genome)))
#print("misc_feature \t" + str(GC_content_misc) + " \t" + str(GC_content_misc/len(genome)))
#print("repeat_region \t" + str(GC_content_repeats) + " \t" + str(GC_content_repeats/len(genome)))
#print("rRNA \t" + str(GC_content_rrna) + " \t" + str(GC_content_rrna/len(genome)))
#print("trna \t" + str(GC_content_trna) + " \t" + str(GC_content_trna/len(genome)))

#"{:%.1f}%.format(GC_content_cds)/(len(genome))")

#to organize this genome..
# dictionaries!
# genes
# key = gene name
# value = [key = exon#, value = seq]

							
	


	
		
		
		
		
	
	
