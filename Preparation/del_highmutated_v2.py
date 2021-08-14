from Bio import AlignIO,SeqIO
import sys

my_fasta = sys.argv[1]
my_out = ".".join(my_fasta.split(".")[:-1])+"_minusbad.fasta"
a=AlignIO.read(my_fasta,"fasta")

cols=[]
for c in range(0,a.get_alignment_length()):
    cols.append(c)

print("indexing columns with unique substitutions...")
col=[]
nt=[]
for c in cols:
    colu=a[:,c]
    As=[]
    Ts=[]
    Gs=[]
    Cs=[]
    Ambs=[]
    for s in colu[::]:
        if s == "A":
            As.append(s)
        elif s == "T":
            Ts.append(s)
        elif s == "G":
            Gs.append(s)
        elif s == "C":
            Cs.append(s)
        else:
            Ambs.append(s)
    if len(As) == 1 and len(Ts) != 1 and len(Gs) != 1 and len(Cs) != 1:
        col.append(c)
        nt.append("A")
    elif len(As) != 1 and len(Ts) == 1 and len(Gs) != 1 and len(Cs) != 1:
        col.append(c)
        nt.append("T")
    elif len(As) != 1 and len(Ts) != 1 and len(Gs) == 1 and len(Cs) != 1:
        col.append(c)
        nt.append("G")
    elif len(As) != 1 and len(Ts) != 1 and len(Gs) != 1 and len(Cs) == 1:
        col.append(c)
        nt.append("C")
    else:
        pass
eval_cols={"col":col,"nt":nt}

print("deleting sequences with more than 0.05 % unique mutations...")
new_fasta=[]
for x,r in enumerate(a):
    print(x)
    pos=[]
    for item,i in enumerate(eval_cols["col"]):
        if r[int(i)] == eval_cols["nt"][item]:
            pos.append(1)
    if len(pos) <= len(cols)*0.0005:
        new_fasta.append(r)
    else:
        pass

print("writing the final file...")
SeqIO.write(new_fasta,my_out,"fasta")
