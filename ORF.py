import File_Path

""" This code below takes Fasta file as an input and find Open Reading Frame (possible protein coding sequences (substrings) of a SEQUENCE (String)) 
    for all records (SEQUENCES). This code considers all three reading frames on a forward strand and give simple analytics for each record, each frame 
    and across all frames providing the details of longest and Shortest ORFs including thier postions on a SEQUENCE, sequence and length."""

def find_ORF(my_fasta):
    from Bio import SeqIO #import SeqIO from Bio package of Python
    seqdict = SeqIO.to_dict(SeqIO.parse(my_fasta, "fasta")) # create a dictionary of sequences with record ID as a key

    frame_with_their_max_lengths = {} #for compiling dictionary of frames with their max Lengths
    frames_of_max_lengths_with_id_positions_sequence = {} #for compiling dictionary of frames with their IDs having Max ORF Seqs and Positions
    
    frame_with_their_min_lengths = {} #for compiling dictionary of frames with their min Lengths
    frames_of_min_lengths_with_id_positions_sequence = {} #for compiling dictionary of frames with their IDs having min ORF Seqs and Positions

    for frame in range(3):  # consider all three reading frames
        ORF_Max_Lengths_with_id = {}
        Longest_ORF_and_Positions_with_id = {}
        
        ORF_Min_Lengths_with_id = {}
        Shortest_ORF_and_Positions_with_id = {}
        
        print(f"*********ANALYZING READING FRAME {frame + 1}***************")

        for record in seqdict:  # iterate through each record in the sequence dictionary
        
            sequence = str(seqdict[record].seq[frame:])  # get the sequence for the current reading frame
            ORF_Positions_and_lenghts_in_a_sequence = {} #for storing positions starting from 0 and useful for slicing the sequence
            

            for nuc in range(0, len(sequence), 3):  # find start codon
                if sequence[nuc:nuc+3] == "ATG":
                    starting_position = nuc
                    for i in range(starting_position, len(sequence), 3):  # scan for stop codon
                        codon = sequence[i:i+3]
                        if codon in ["TAG", "TAA", "TGA"]:  
                            end_position = i + 3
                            ORF_length = end_position - starting_position
                            
                            ORF_Positions_and_lenghts_in_a_sequence[(starting_position, end_position)] = ORF_length
                            
                            break  # stop scanning after finding the stop codon for this ORF in asequence
        
            #outputting results for particular record            
            print("\n")
            print(f"----------> ORF ANALYSIS FOR RECORD ID: {record} IN READING FRAME {frame + 1}")
            print("\n")

            
            Sequence_and_Position_max = []
            Sequence_and_Position_min = []                
            for positions, length in ORF_Positions_and_lenghts_in_a_sequence.items():

                if length == max(ORF_Positions_and_lenghts_in_a_sequence.values()):  # for longest ORF in a record
                                        
                    sequence_with_max_ORF_length = (sequence[positions[0]:positions[1]])
                    Actual_positions_of_longest_ORF = (positions[0] + (frame + 1), positions[1] + frame) #Extract Actual Positions of ORF from an original sequence with index starting at 1     
                    Sequence_and_Actual_Positions_max_length = ("Sequence:", sequence_with_max_ORF_length, "position:", Actual_positions_of_longest_ORF)
                    Sequence_and_Position_max.append(Sequence_and_Actual_Positions_max_length) # for appending sequences and their positions when more than one sequences with maximum length, if any.

                    ORF_Max_Lengths_with_id[record] = length #Compile a dictionary of max ORF lengths with record ids
                    Longest_ORF_and_Positions_with_id[record] = Sequence_and_Position_max #Compile a dictionary of longest ORFs and thier positions with record ids
                    print("Longest ORF and actual Position:", Sequence_and_Actual_Positions_max_length, "Length:", length)
                    print("\n")

                if length == min(ORF_Positions_and_lenghts_in_a_sequence.values()): # for shortest ORF in a record
                    
                    sequence_with_min_ORF_length = (sequence[positions[0]:positions[1]])
                    Actual_positions_of_shortest_ORF = (positions[0] + (frame + 1), positions[1] + frame) #Extract Actual Positions of ORF from an original sequence with index starting at 1
                    Sequence_and_Actual_Positions_min_length = ("Sequence:", sequence_with_min_ORF_length, "position:", Actual_positions_of_shortest_ORF)
                    Sequence_and_Position_min.append(Sequence_and_Actual_Positions_min_length) ## for appending sequences and their positions when more than one sequences with minumum length, if any

                    ORF_Min_Lengths_with_id[record] = length #Compile a dictionary of min ORF lengths with record ids
                    Shortest_ORF_and_Positions_with_id[record] = Sequence_and_Position_min #Compile a dictionary of shortest ORFs and their positions with record ids
                    print("Shortest ORF and actual Position:", Sequence_and_Actual_Positions_min_length, "Length:", length)
                    print("\n")        

        #Outputting the results for particular Frame

        print(f"*************************************\n##################.....SUMMARY OF ORF ANALYSIS ACROSS ALL RECORDS in FRAME {frame + 1}.....##################\n*************************************")
        print("\n")
        print(f"Record ID (s) having maximum ORF length in Frame :{frame + 1}:")

        records_with_position_Seq_Max = []  # A list to store record IDs having longest ORFs with thier sequence and positions 
        records_with_position_Seq_Min = []  # A list to store record IDs having shortest ORFs with thier sequence and positions

        for record_id, max_length in ORF_Max_Lengths_with_id.items():
            if max_length == max(ORF_Max_Lengths_with_id.values()):

                record_and_thier_seq_pos_max = (record_id, Longest_ORF_and_Positions_with_id[record_id]) #club record IDs with their longest ORFs and their Positions
                records_with_position_Seq_Max.append(record_and_thier_seq_pos_max) # Append Record IDs when more than one record IDs having longest ORFs
                
                frame_with_their_max_lengths[frame + 1] = max_length
                frames_of_max_lengths_with_id_positions_sequence[frame + 1] = records_with_position_Seq_Max

                print(f"{record_id}, Longest ORF Length: {max_length}, Longest ORF Sequence and Actual Position: {Longest_ORF_and_Positions_with_id[record_id]}")                
                print("\n")

        print("\n")

        print(f"Record ID (s) having minimum ORF length in Frame: {frame + 1}")

        for record_id, min_length in ORF_Min_Lengths_with_id.items():
            if min_length == min(ORF_Min_Lengths_with_id.values()):
                record_and_their_seq_pos_min = (record_id, Shortest_ORF_and_Positions_with_id[record_id]) #club record IDs with their shortest ORFs and their Positions
                records_with_position_Seq_Min.append(record_and_their_seq_pos_min) # Append Record IDs when more than one record IDs having shortest ORFs
                
                frame_with_their_min_lengths[frame + 1] = min_length
                frames_of_min_lengths_with_id_positions_sequence[frame + 1] = records_with_position_Seq_Min
                print(f"{record_id}, Shortest ORF Length: {min_length}, Shortest ORF Sequence and Actual Position: {Shortest_ORF_and_Positions_with_id[record_id]}")                
                print("\n")


    print("################################################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n##################.....SUMMARY OF ORF ANALYSIS ACROSS ALL FRAMES.....################################################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n##################")
    print("\n")
    print("\n")

    # Outputting the results for all the frames with their summary

    print("Frames, IDs, Position and sequence having maximum length across all records and across all frames")
    for frame, max_length in frame_with_their_max_lengths.items():
        if max_length == max(frame_with_their_max_lengths.values()):
            print(f"Frame:{frame},\n maximum length:{max_length}\n, IDs, Sequence and Positions: {frames_of_max_lengths_with_id_positions_sequence[frame]}")
    print("\n")

    print("Frames, IDs, Position and sequence having minimum length across all records and across all frames")
    for frame, min_length in frame_with_their_min_lengths.items():
        if min_length == min(frame_with_their_min_lengths.values()):
            print(f"Frame:{frame},\n minimum length:{min_length}\n, IDs, Sequence and Positions: {frames_of_min_lengths_with_id_positions_sequence[frame]}")



    print("*************************************\n##################.....END OF ORF ANALYSIS.....##################\n*************************************")





if __name__ == "__main__": #Picks fasta file of user selection and run the mentioned function on it.
    fasta_path = File_Path.pick_fasta_file()
    if fasta_path:
        find_ORF(fasta_path)
    else:
        print("No file selected.")
        