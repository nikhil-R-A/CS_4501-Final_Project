#!usr/bin/env python

import os
import sys
import re
from Bio import SeqIO

my_fasta = sys.argv[1]
p1 = int(sys.argv[2])+1

def extend_str(s, l):
    return (s*l)[:l]

my_out = my_fasta[:-6]+"_"+str(p1-1)+"Ns.fasta"
c_number_of_Ns_removed = 0
my_filter_r = []
for r in SeqIO.parse(my_fasta,"fasta"):
	my_seq = str(r.seq)
	Ns_to_check = extend_str("N",p1)
	my_re_comp = re.compile(re.escape(Ns_to_check), re.IGNORECASE)
	my_re_find = re.findall(my_re_comp,my_seq)
	#print(Ns_to_check)
	#print(my_re_find)

	if len(my_re_find) == 0:
		my_filter_r.append(r)
		#print("selected",my_seq)
	else:
		#print("removed",my_seq)
		c_number_of_Ns_removed += 1

SeqIO.write(my_filter_r,my_out,"fasta")
print("Sequences removed: ",c_number_of_Ns_removed)
print("Sequences selected in ",my_out," :",len(my_filter_r))
	