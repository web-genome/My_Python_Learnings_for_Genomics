""" This code below returns the reulsts of BLAST over all the sequences of FASTA File. """

import File_Path

def test_blast(my_fasta):
    import Bio
    from Bio import Blast
    from Bio.Blast import NCBIWWW, NCBIXML

    from Bio import SeqIO #import SeqIO from Bio package of Python
    seqdict = SeqIO.to_dict(SeqIO.parse(my_fasta, "fasta")) # create a dictionary of sequences with record ID as a key

    for record in seqdict:

        sequence = str(seqdict[record].seq)

        result_handle = NCBIWWW.qblast("blastn", "nt", sequence) #Run qblast on sequence of particular record by searching for nucleotide database
        blast_record = NCBIXML.read(result_handle) #Read and store results in XML format
        for alignment in blast_record.alignments: #check for alignments in blast records
            for hsp in alignment.hsps:
                if hsp.expect < 0.01: #check for high scoring pair with less than 1 % score
                    print("****Alignment****")
                    print (f"*******{record}*********")
                    print("sequence:", alignment.title)
                    print("length:", alignment.length)
                    print("e value:", hsp.expect)
                    print(hsp.query)
                    print(hsp.match)
                    print(hsp.sbjct)
                    print("\n")


    return test_blast(my_fasta)

if __name__ == "__main__": #Picks fasta file of user selection and run the mentioned function on it.
    fasta_path = File_Path.pick_fasta_file()
    if fasta_path:
        test_blast(fasta_path)
    else:
        print("No file selected.")