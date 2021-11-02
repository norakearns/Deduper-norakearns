#!/usr/bin/env python


import argparse
import sys
import re

def get_args():
    parser = argparse.ArgumentParser(description="A program to eliminate PCR duplicates from a sam file")
    parser.add_argument("-f", help="input sam file of aligned reads", required=True, type=str)
    parser.add_argument("-p", help="designates file is paired-end", required=False, action='store_true')
    parser.add_argument("-o", help="output file with deduplicated reads", required=False, default="sam_deduped.sam", type=str)
    parser.add_argument("-u", help="umi file with list of known umis", required=False, default="STL96.txt", type=str)
    return parser.parse_args()
     

args = get_args()

Sam_f = args.f
Out_sam = args.o
Umi_f = args.u

if args.p: # throw an error if paired-end option used, and exit the program
    print("Error: Paired-end not implemented for this program")
    sys.exit(1)

Sam_file = open(Sam_f, "r") # opens the argument (-f, non-deduplicated input sam file) for reading
Umi_file = open(Umi_f, "r") # opens the argument (-u, umi file) for reading
Output_sam = open(Out_sam, "wt") # opens the argument (-o, output file) for writing
Wrong_umi = open("wrong_umi_file.sam", "wt") # opens a file to write the reads with incorrect errors for writing

Umi_list = []
for line in Umi_file:
    # make a list of known UMIs
    stripped_line = line.rstrip('\n')
    Umi_list.append(stripped_line) 

def FWD_Correct_pos(CigarString, Pos):
    '''
    Takes as it's input a cigar string and the left-most mapping position, and edits the position to be the true, unclipped start position 
    For a forward read: take the number of nucleotides the read is soft-clipped by and subtract it from the position.
    '''

    if "S" in CigarString[:-1]: # Only get 5' soft clipping
        Num_Clipped = CigarString.split("S")[0] # get the number of soft clipped reads by splitting before the "S." Ex: (2, S71M)
        True_Pos = int(Pos) - int(Num_Clipped) # subtract soft-clipped nt from the position
    else:
        True_Pos = Pos
    return str(True_Pos)

def REV_Correct_pos(CigarString, Pos, Read_len):
    '''
    For a reverse read: 
       Position
       -  # of nucleotides the read is soft-clipped by at the left (3') end (S)
       -  # of nucleotides inserted (I)
       +  # of nucleotides deleted (D)
       +  # of nucleotides gapped (N)
       + Read length
       -------------
       True Position
       '''
    cigar_dict = {} # dictionary to hold items in cigar string. Structure ex: {'S':[2], 'I':[3,4], 'D':[2],'N':[1024]}
    pattern = '([0-9]+)([IND])' # make a pattern for each insertion, deletion, or gap occurrence in the cigar string. Ex: 2N, 3I
    s_pattern = '([0-9]+)([S])' # make a pattern for the soft-clipping. Ex: 2S
    matches = re.findall(pattern, CigarString) # use regular expression matching to find all occurences of insertions, deletions, and gaps
    s_matches = re.findall(s_pattern, CigarString[:-1]) # use regex to find soft clipping at the 3' end
    for item in matches: # Ex: [(2,I),(3,I)]
        if item[1] in cigar_dict.keys(): # If one occurrence is already in there
            cigar_dict[item[1]].append(int(item[0])) # Append the number from the second occurrence
        else:
            cigar_dict[item[1]]=[int(item[0])] # Add a k,v pair to the dictionary with the form {'I':[2]}
    for s_item in s_matches:
        cigar_dict[s_item[1]]=[int(s_item[0])] # Add the 3' soft clipping to the dictionary in the form {'S':[2]}
    if 'S' in cigar_dict.keys():
        Num_Clipped = cigar_dict['S'][0] # Grab the value from the dictionary for the key 'S'
    else:
        Num_Clipped = 0 # if no occurrences of S
    if 'I' in cigar_dict.keys():
        Num_Inserted = sum(cigar_dict['I']) # sum insertions
    else:
        Num_Inserted = 0
    if 'N' in cigar_dict.keys():
        Num_Gapped = sum(cigar_dict['N']) # sum gaps
    else:
        Num_Gapped = 0
    if 'D' in cigar_dict.keys():
        Num_Deleted = sum(cigar_dict['D']) # sum deletions
    else:
        Num_Deleted = 0

    adjust = Read_len - Num_Clipped - Num_Inserted # what to adjust the position by
    True_pos = int(Pos) + adjust + Num_Gapped + Num_Deleted
    return str(True_pos)


Read_Dict = {} # Make a dictionary of the form {'UMI-Chrom-pos-strand': sam line}, this dict will be cleared out after each chromosome
last_chrom = '0'
for line in Sam_file:
    Read_array = [] # Make an array with all of the important identifiers from each read [UMI, chrom, pos, strand]
    if line.startswith("@"):
        Output_sam.write(line)
        continue
    entry = line.split('\t') # split the line on tabs to make an list of each line
    chrom = entry[2] # get chrom from position 3 in the array
    cigar = entry[5] # get cigar string from position 6 in the array
    Read_len = len(entry[9]) # get the read length 
    # if chrom is not equal to the previous one, write out the contents of the dictionary to keep the dictionary small and reduce lookup time
    if chrom != last_chrom: 
        last_chrom = chrom # replace the "previous" chrom value with the new one
        for value in Read_Dict.values(): # write out the keys in the dictionary (the sam lines) to the output file
            Output_sam.write('{}'.format(value))
        Read_Dict = {} # empty out the dictionary
    qname_umi = entry[0][-8:] # grab the umi from the Qname
    if qname_umi in Umi_list: 
        Read_array.append(qname_umi)
    else:
        Wrong_umi.write(line) # if the umi isn't in the known list, write the corresponding sam line out to a "wrong umi" file
        continue
    Read_array.append(chrom) 
    pos = entry[3] # get the mapping position from position 4 in the list
    flag = int(entry[1])
    if ((flag & 16) == 16):  # if the bit at 16 is set, read is on the reverse strand
        Read_array.append("Rev")
        True_Pos = REV_Correct_pos(cigar, pos, Read_len) # Correct the mapping position
    else:
        Read_array.append("Fwd")
        True_Pos = FWD_Correct_pos(cigar, pos) # Correct the mapping position
    Read_array.append(True_Pos) 
    Read_string = Read_array[0] + "_" + Read_array[1] + "_" + Read_array[2] + "_" + Read_array[3] # Concatenate the items in the array to a string: 'UMI-chrom-pos-strand'
    if Read_string in Read_Dict: # if the 'UMI-chrom-pos-strand' string already exists as a key in the dictionary, don't add the current read
        continue
    else:
        Read_Dict[Read_string] = line # if the 'UMI-chrom-pos-strand' string doesn't exist as a key in the dictionary, add the current read

for value in Read_Dict.values(): # when you get through the last entry of the last chromosome, empty the contents of the dictionary 
    Output_sam.write('{}'.format(value))