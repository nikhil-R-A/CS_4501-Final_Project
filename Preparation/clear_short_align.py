import sys
from Bio import SeqIO

fasta_file=sys.argv[1]
min_length=sys.argv[2]
my_out = ".".join(fasta_file.split(".")[:-1])+"_align_clear.fasta"

sequences=[]
for i in SeqIO.parse(fasta_file, "fasta"):
    count=0
    for x in i.seq:
        if x == "-":
            count=count+1
    if float(count) <= float(len(i.seq)*float(min_length)):
        sequences.append(i)
             
SeqIO.write(sequences,my_out,"fasta")
print("Finish!")  