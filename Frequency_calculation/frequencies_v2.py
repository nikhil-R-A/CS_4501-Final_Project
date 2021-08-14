import time
import sys
start_time = time.time()

my_aln = sys.argv[1]
my_aln_open = open(my_aln,"r")

# retrieving the sequences of the fasta file
my_seq = ""
my_sequences = []
for line in my_aln_open:
    line = line.strip()    
    start_seq = False
    if line.startswith(">"):
        my_sequences.append(my_seq)
        my_seq = ""
    else:
        start_seq = True

    while start_seq:
        my_seq += line
        start_seq = False
my_sequences_all = (my_sequences[1:]) + [my_seq]

print("number of sequences: ",len(my_sequences_all))

# List of hashes for each colum
hl = []
for i in range(len(my_sequences_all[0])):
    hl.append({'A': 0, 'C': 0, 'G': 0, 'T': 0, '-' : 0, 'N': 0, 'M': 0, 'R': 0, 'W': 0, 'S' : 0,\
               'Y': 0, 'K': 0, 'V': 0,  'H': 0, 'D': 0, 'B': 0})

# Couting the bases
nLines = 0
for line in my_sequences_all:
    for idx, c in enumerate(line):
        hl[idx][c] += 1
    nLines += 1

results = []
nLines = float(nLines)
for char in ['A', 'C', 'G', 'T' , '-', 'N', 'M', 'R', 'W' , 'S', \
             'Y', 'K', 'V', 'H' , 'D', 'B']:
    freq =  [round(x[char]/nLines,20) for x in hl]
    #row = str("{}\t{}".format(char, "\t".join(["{:0.2f}".format(x[char]/nLines) for x in hl])))
    if len(set(freq)) == 1 and list(set(freq))[0] == 0.0:
        pass
    else:
        results.append([char]+list(map(str,freq)))

# transpoing list of list
results_tansp = []
for i in zip(*results):
    results_tansp.append(list(i))

count_positio = 0
results_tansp2 = [] 
header = ["Position"] + results_tansp[0]
for ii in results_tansp[1:]:
    results_tansp2.append([str(count_positio)]+ii)
    count_positio += 1
results_tansp2.insert(0,header)

# writing the results
my_results_open = open(my_aln+".freq","w",)
for r in results_tansp2:
    my_results_open.write("\t".join(r))
    my_results_open.write('\n')
    
my_aln_open.close()
my_results_open.close()
print("End ... check: ",my_aln+"_rel.freq")
print("--- %s seconds ---" % (time.time() - start_time))