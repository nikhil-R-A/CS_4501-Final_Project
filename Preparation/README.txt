1.- Enter to GISAID (gisaid.org) and download the sequences and the metadata.
	1.1.- Login
	1.2.- Click en la opción "downloads"
	1.3.- En la sección Download packages descargar, FASTA y metadata
2.- A File named "sequences_fasta_2021_05_25.tar.xz" and "metadata_tsv_2021_05_25.tar.xz" will be obtained
3.- Decompress the two files to obtain the files: "metadata.tsv" and "sequences.fasta"
tar -xf sequences_fasta_2021_05_25.tar.xz && tar -xf metadata_tsv_2021_05_25.tar.xz
4.- Remove sequences with less than 29 000 nt using the script "clear_short_seqs.py".
	4.1.- Execute the next command line: python3 clear_short_seqs.py sequences.fasta 29000
	4.2.- output=sequences_clear.fasta
5.- Align the sequences using ViralMSA.py using genome "EPI_ISL_406801.fa" as the reference:
ViralMSA.py -s sequences_clear.fasta -r EPI_ISL_406801.fa -e santiago.jus.are@usp.br -o alignment -t 20 --omit_ref
	5.6.- A folder named alignment will be created with the file named: sequences_clear.fasta.aln, change the name to "aligned.fasta"
6.- Transfer the file "aligned.fasta" to the previous folder 
7.- Remove sequences with more than 290 Ns using the script "remove_Ns.py" 
command line: python3 remove_Ns.py aligned.fasta 290
output= aligned_290Ns.fasta
8.- Remove sequences with more than 0.05 % unique mutations using the script "del_highmutated_v2.py"
command line: python3 del_highmutated_v2.py aligned_290Ns.fasta
output= aligned_290Ns_minusbad.fasta
9.- Remove sequences with more than 2 % gaps using the script "clear_short_align.py"
command line: python3 clear_short_align.py aligned_290Ns_minusbad.fasta 0.02
output= aligned_290Ns_minusbad_align_clear.fasta