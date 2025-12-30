#This test computes all the below functions on given Fasta File simultaneously.

import File_Path
from Record_length import record_length
from sequence_length import seq_length
from ORF import find_ORF
from kmer import find_kmers
from Blast_alignment import test_blast
import sys

if __name__ == "__main__":
    fasta_path = File_Path.pick_fasta_file()
    if fasta_path:
        
        sys.stdout.flush()
        find_kmers(fasta_path) #This function call for user defined k-mer length 
        sys.stdout.flush()

        
        sys.stdout.flush()
        seq_length(fasta_path)
        sys.stdout.flush()

        
        sys.stdout.flush()
        find_ORF(fasta_path)
        sys.stdout.flush()

        
        sys.stdout.flush()
        record_length(fasta_path)
        sys.stdout.flush()

        sys.stdout.flush()
        test_blast(fasta_path)
        sys.stdout.flush()

    else:
        print("No file selected.")
        sys.stdout.flush()