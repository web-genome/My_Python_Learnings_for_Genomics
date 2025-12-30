
import File_Path   
"""This code returns total records contained in the Fatsa file."""

def record_length(my_fasta):
    from Bio import SeqIO #import SeqIO from Bio package of Python
    records = list(SeqIO.parse(my_fasta, "fasta")) # create a dictionary of sequences with record ID as a key
    print("*********TOTAL RECORDS***************")
    print("\n")
    print("There are %i records found in the input file" % len(records)) # gets the number of records in the file as they are identified using the keys represented by each record description

if __name__ == "__main__": #Picks fasta file of user selection and run the mentioned function on it.
    fasta_path = File_Path.pick_fasta_file()
    if fasta_path:
        record_length(fasta_path)
    else:
        print("No file selected.")
# Record length.py


