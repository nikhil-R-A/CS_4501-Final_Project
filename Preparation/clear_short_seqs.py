import sys
from Bio import SeqIO

fasta_file=sys.argv[1]
min_length=sys.argv[2]
my_out = ".".join(fasta_file.split(".")[:-1])+"_clear.fasta"

sequences=[]
for i in SeqIO.parse(fasta_file, "fasta"):
    if (len(i.seq) >= float(min_length)):
        sequences.append(i)

             
SeqIO.write(sequences,my_out,"fasta")
print("Finish!")  