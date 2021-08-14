#!usr/bin/env python

import os
import sys
import csv
import collections
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

my_aln = sys.argv[1]#".aln"
my_pos = sys.argv[2]# .txt

dic_positions = {}
for pos in csv.reader(open(my_pos),delimiter="\t"):
	dic_positions[pos[0]] = pos[1]
print(dic_positions)
recs_alin  = []
size2check = []
for rec in SeqIO.parse(my_aln,"fasta"):
	recs_alin.append(rec)
	size2check.append(len(rec.seq))

recs_seq_selec = []
recs_seq_selec_2 = []
id_list = []
if len(set(size2check)) == 1:
	for r in recs_alin:
		my_flags = []
		to_check = []
		for k,v in dic_positions.items():
			if r[int(k)-1] == str(v).lower() or r[int(k)-1] == str(v).upper():
				my_flags.append("True")
				to_check.append(str(int(k)-1)+":"+r[int(k)-1]+"(True)")
			else:
				my_flags.append("False")
				to_check.append(str(int(k)-1)+":"+r[int(k)-1]+"(False)")
		if len(set(my_flags)) == 1 and list(set(my_flags))[0] == "True":
			id_list.append(r.id)
			recs_seq_selec_2.append(r)
			rec_2 = SeqRecord(Seq(str(r.seq).replace("-",""),IUPAC.IUPACUnambiguousDNA),id=r.id,description=r.description)
			recs_seq_selec.append(rec_2)
SeqIO.write(recs_seq_selec,my_aln[:-6]+"_"+my_pos[:-4]+"_seq.fasta","fasta")
SeqIO.write(recs_seq_selec_2,my_aln[:-6]+"_"+my_pos[:-4]+"_aligned.fasta","fasta")
with open(my_aln[:-6]+"_"+my_pos[:-4]+".list","w") as my_list:
	w_list = csv.writer(my_list,delimiter=",",quotechar='\"', quoting=csv.QUOTE_NONNUMERIC)
	w_list.writerow(id_list)