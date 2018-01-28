import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref

allowedMismatch = 5
hashLength = 2
def makeHashTable(ref):
	#reference = read_reference(reference_fn)

	hashTable = {}
	i = 0
	while i < len(ref) - hashLength + 1:
		try:
			hashTable[ref[i:i+hashLength]].append(i)
		except:
			hashTable[ref[i:i+hashLength]] = []
			hashTable[ref[i:i+hashLength]].append(i)
		i = i + 1
	return hashTable



def betterAlgo(paired_end_reads, ref, hashTable):
  	all_read_alignment_locations = []
  	output_readPairs = []
  	count = 0
  	start = time.clock()
  	for read_pair in paired_end_reads:
  		print "read_pair " + read_pair
  		count += 1
  		read_alignment_locations = []
  		output_read_pair = []
  		if count % 10 == 0:
  			time_passed = (time.clock()-start)/60
  			print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
  			remaining_time = time_passed/count*(len(paired_end_reads)-count)
  			print 'Approximately {:.3} minutes remaining'.format(remaining_time)
  		for read in read_pair:
  			print "for read " + read
  			min_mismatches = len(read) + 1
  			min_mismatch_location = -1
           	#check regular version of this read
           	print len(read)/hashLength
           	for i in range(len(read)/hashLength):
           		print "for i in range"
           		startIndexRead = (i*hashLength)
           		hashSnip = read[startIndexRead:((i + 1)*hashLength)]
           		try:
           			possibleIndexList = hashTable[hashSnip]
           		except: 
           			print "no match for " + hashSnip
           			continue
           		for possIndex in possibleIndexList:
           			print "possIndex " + str(possIndex)
           			if possIndex >= startIndexRead:
           				#not out of range
           				startReadCompare = 0
           				startRefCompare = possIndex - startIndexRead
           				print "snip " + read
           				print "ref snip " + ref[startReadCompare:startReadCompare + len(read)]
           				n_mismatches = 0
           				for i in range(len(read)):
           					if read[i] != ref[startRefCompare + i]:
           						n_mismatches = n_mismatches + 1

           	#then check reversed


	return all_read_alignment_locations, output_readPairs

reference = "AATTCCGGAATTCCGG"
#reference = read_reference(reference_fn)
test = makeHashTable(reference)
reads = ["AAG", "CCT"]
betterAlgo(reads, reference, test)
print test

# if __name__ == "__main__":
#     data_folder = 'practice_W_1'
#     input_folder = join('../data/', data_folder)
#     f_base = '{}_chr_1'.format(data_folder)
#     reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
#     start = time.clock()
#     input_reads = read_reads(reads_fn)
#     # This will take a while; you can use an array slice for example:
#     #
#     #   input_reads = reads[:300]
#     #
#     # to generate some data quickly.

#     reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
#     reference = read_reference(reference_fn)
#     alignments, reads = betterAlgo(input_reads, reference)#trivial_algorithm(input_reads, reference)
#     print alignments
#     print reads
#     output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
#     output_fn = join(input_folder, 'aligned_{}.txt'.format(f_base))
#     with(open(output_fn, 'w')) as output_file:
#         output_file.write(output_str)




