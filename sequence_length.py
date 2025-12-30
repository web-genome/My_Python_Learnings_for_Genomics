import File_Path   

"""This code provides details of genomic sequences contained in the Fasta File pertaining to their lengths."""

def seq_length(my_fasta):
    import pprint
    from Bio import SeqIO #import SeqIO from Bio package of Python
    seqdict = SeqIO.to_dict(SeqIO.parse(my_fasta, "fasta")) # create a dictionary of sequences with record ID as a key
    
    Sequence_Lengths_with_id = {}
    for record in seqdict:  # iterate through each record in the sequence dictionary
        Sequence_Lengths_with_id[record] = len(seqdict[record].seq)  # get length of each sequence and its corresponding id
    print("************************LENGTHS WISE DISPLAY************************")
    print("\n")
    print("Sequence Lengths with IDs:")
    pprint.pprint(Sequence_Lengths_with_id)
    print("\n")

    # Calculate min and max lengths of all sequences and their counts
    min_length = min(Sequence_Lengths_with_id.values())
    max_length = max(Sequence_Lengths_with_id.values())
    total_number_of_minimum_length_sequences = list(Sequence_Lengths_with_id.values()).count(min_length)
    total_number_of_maximum_length_sequences = list(Sequence_Lengths_with_id.values()).count(max_length)
    print(f"Total number of sequences with minimum length: {total_number_of_minimum_length_sequences}")
    print(f"Total number of sequences with maximum length: {total_number_of_maximum_length_sequences}")
    print("\n")

    for seq_id, length in Sequence_Lengths_with_id.items():
        if length == min_length:
            print(f"ID_of_seq_with_min_length: {seq_id}, Length: {length}") # print the id, length and sequence of the shortest sequences
            print(f"Sequence: {seqdict[seq_id].seq}")
        else:
            continue
    print("\n")

    for seq_id, length in Sequence_Lengths_with_id.items():        
        if length == max_length:
            print(f"ID_of_seq_with_max_length: {seq_id}, Length: {length}") # print the id, length and sequence of the longest sequences
            print(f"Sequence: {seqdict[seq_id].seq}")
        else:
            continue



if __name__ == "__main__": #Picks fasta file of user selection and run the mentioned function on it.
    fasta_path = File_Path.pick_fasta_file()
    if fasta_path:
        seq_length(fasta_path)
    else:
        print("No file selected.")


