import File_Path

""" This code takes Fasta File as an input and ask user to provide an integer "k" as another input. This integer indicates the
    length of the k-mer (substring of any genomic sequence) to look for and finds the most occuring k-mer of length "k" in 
    for each sequence and across all the sequences."""

def find_kmers(my_fasta):

    from Bio import SeqIO #import SeqIO from Bio package of Python
    seqdict = SeqIO.to_dict(SeqIO.parse(my_fasta, "fasta")) # create a dictionary of sequences with record ID as a key
    print("\nAll previous outputs completed. Proceeding to k-mer length input.")

    import sys
    sys.stdout.flush()

    kmer_length = int(input("Enter the length of the k-mer: ")) # define the length of the k-mer from a user input
   
    kmer_counts_in_all_sequences = {} # dictionary to store k-mer counts across all sequences

    for record in seqdict:  # iterate through each record in the sequence dictionary
        sequence = seqdict[record].seq        
        kmer_counts_in_sequence = {}  # dictionary to store k-mer counts of a sequence contained in a record
        
        for i in range(len(sequence) - kmer_length + 1): # iterate through the sequence to extract k-mers and count their occurrences
            kmer = str(sequence[i:i + kmer_length])
            kmer_counts_in_sequence[kmer] = 0
        
        for i in range(len(sequence) - kmer_length + 1):
            kmer = str(sequence[i:i + kmer_length])
            if kmer == kmer:
                kmer_counts_in_sequence[kmer] += 1
                               

        for kmer, count in kmer_counts_in_sequence.items(): # iterate through the k-mer counts to find the most frequent k-mer in a sequence
            
            if count == max(kmer_counts_in_sequence.values()):
                print("\n")
                print(f"---------------K-MER ANALYSIS FOR RECORD ID: {record} ----------------")
                print("\n")
                print(f"Most frequent k-mer in record ID {record}: {kmer}, Count: {count}")
                print("\n")

            if kmer in kmer_counts_in_all_sequences: #store each kmer and their counts from each sequence to the dictionary for all sequences 
                kmer_counts_in_all_sequences[kmer] += count
            else: 
                kmer_counts_in_all_sequences[kmer] = count
        

            
    print("*********SUMMARY OF K-MER COUNTS ACROSS ALL RECORDS***************")
    print("\n")
    for kmer, total_count in kmer_counts_in_all_sequences.items():
        if total_count == max(kmer_counts_in_all_sequences.values()): # find the most frequent k-mer across all sequences
            print(f"Most frequent k-mer across all records: {kmer}, Total Count: {total_count}")



if __name__ == "__main__": #Picks fasta file of user selection and run the mentioned function on it.
    fasta_path = File_Path.pick_fasta_file()
    if fasta_path:
        find_kmers(fasta_path)
    else:
        print("No file selected.")   


