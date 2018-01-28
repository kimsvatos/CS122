import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref

hashLength = 10
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
  	output_read_pairs = []
  	count = 0
  	start = time.clock()
  	for read_pair in paired_end_reads:
  		#print "read_pair " + read_pair
  		count += 1
  		read_alignment_locations = []
  		output_read_pair = []
  		if count % 10 == 0:
  			time_passed = (time.clock()-start)/60
  			#print '{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed)
  			remaining_time = time_passed/count*(len(paired_end_reads)-count)
  			#print 'Approximately {:.3} minutes remaining'.format(remaining_time)
  		for read in read_pair:
  			#print "for read in readpair" + read
  			#rangeNum = len(read) / hashLength
  			
  			min_mismatches = len(read) + 1
  			min_mismatch_location = -1
  			for i in range(len(read)/hashLength):
  			#while i < len(read)/hashLength:
				startIndexRead = (i*hashLength)
				hashSnip = read[startIndexRead:((i + 1)*hashLength)]
				#print "hashSnip " + hashSnip
				try:
					possibleIndexList = hashTable[hashSnip]
				except: 
					#print "no match for " + hashSnip
					continue
				for possIndex in possibleIndexList:
					#print "possIndex " + str(possIndex)
					if possIndex >= startIndexRead and (possIndex - startIndexRead + len(read)) < len(ref):
           				#not out of range
						startReadCompare = 0
						startRefCompare = possIndex - startIndexRead
						#print "snip " + read
						#print "ref snip " + ref[startRefCompare:startRefCompare + len(read)]
						n_mismatches = 0
						for i in range(len(read)):
							print "len read " + str(len(read))
							print "len ref " + str(len(ref))
							print "read index " + str(i)
							print "index " + str(startRefCompare + i)
							if read[i] != ref[startRefCompare + i]:
								n_mismatches = n_mismatches + 1
						#print "mismatches " + str(n_mismatches)
						if n_mismatches < min_mismatches:
							min_mismatches = n_mismatches
							min_mismatch_location = startRefCompare
							#print "min location: " + str(startRefCompare)

			reversed_read = read[::-1]
			for i in range(len(read)/hashLength):
  			#while i < len(read)/hashLength:
				startIndexRead = (i*hashLength)
				hashSnip = reversed_read[startIndexRead:((i + 1)*hashLength)]
				#print "hashSnip " + hashSnip
				try:
					possibleIndexList = hashTable[hashSnip]
				except: 
					#print "no match for " + hashSnip
					continue
				for possIndex in possibleIndexList:
					#print "possIndex " + str(possIndex)
					if possIndex >= startIndexRead and (possIndex - startIndexRead + len(read)) < len(ref):
           				#not out of range
						startReadCompare = 0
						startRefCompare = possIndex - startIndexRead
						#print "snip " + reversed_read
						#print "ref snip " + ref[startRefCompare:startRefCompare + len(read)]
						n_mismatches = 0
						for i in range(len(read)):
							if reversed_read[i] != ref[startRefCompare + i]:
								n_mismatches = n_mismatches + 1
						#print "mismatches " + str(n_mismatches)
						if n_mismatches < min_mismatches:
							min_mismatches = n_mismatches
							min_mismatch_location = startRefCompare
							read = reversed_read
							#print "min location: " + str(startRefCompare)

			read_alignment_locations.append(min_mismatch_location)
			output_read_pair.append(read)
		all_read_alignment_locations.append(read_alignment_locations)
		output_read_pairs.append(output_read_pair)
           
			
	return all_read_alignment_locations, output_read_pairs

#reference = "AATTCCGGAATTCCGG"
#reference = read_reference(reference_fn)
#test = makeHashTable(reference)
#reads = [["AAGT", "CCTG"], ["CCAA", "CCTG"]]
#a, b = betterAlgo(reads, reference, test)

#print a
#print b

if __name__ == "__main__":
    data_folder = 'practice_W_1'
    input_folder = join('../data/', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)
    #print "PRINT BELOW"
    #test_reads = []
    #for i in range(10):
       # print input_reads[i]
     #   test_reads.append(input_reads[i])
    # This will take a while; you can use an array slice for example:
    #
    #   input_reads = reads[:300]
    #
    # to generate some data quickly.

    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
 	
    alignments, reads = betterAlgo(input_reads, reference, makeHashTable(reference))
    print alignments
    #print reads
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




