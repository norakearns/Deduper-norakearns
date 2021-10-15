
# Assignment Outline: PCR Deduplication
Write up a strategy for writing a Reference Based PCR Duplicate Removal tool. That is, given a sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory. You should not write any code for this portion of the assignment. 

# Define the Problem
PCR is used to amplify many copies of a read prior to sequencing. This is done to improve sequencing quality, as sequencing errors are unlikely to occur at the exact same nucleotide across all duplicates, and the highest quality/consensus nucleotide across all duplicates can be taken as the "true" base call.

However, PCR duplicates must be removed for downstream processing (genome assembly, variant calling, etc.). If they are not removed, they would be treated as "molecular duplicates," or identical sequences which naturally occur in the genome being studied. If assembling a genome, this would result in the same sequence being inserted multiple times, even if that does not truly happen in the genome.  

If PCR duplicates are not removed for alignment to a reference, they can also falsely increase the coverage of a specific region (meaning that instead of having many independent reads overlapping the same reference segment but with different start/end positions, there are multiples of the same exact read overlapping the exact same segment of a reference with the same start/end position).

# What does a PCR-duplicate look like
> UMI 
> Strand (strand specific?) (FLAG, SAM col 2)
> Chromosome (RNAME, SAM col 3)
> Alignment position (POS, SAM col 4): Be aware of Soft clipping (in the CIGAR string)

# Functions
Correct_pos(phred_score: str) -> int
```Takes as it's input a cigar string and the position, and edits the position to be the true, unclipped start position (take the number of nucleotides the read is soft-clipped by and add it to the position).```
    Input: 2S69M
    Input: 76814286
    Return(Corrected Position): 76814284

    * Use Bi621 PS8: ParseSam.py to identify if something is on the forward or reverse strand from the bitwise flag.
        - if ((flag & 4) == 0) # if the 4th bit in the flag is 0, the read is mapped
        - example: 80 = reverse strand
                   143 = forward strand

Create_String(quality_line: str) -> float
```Takes as it's input an whole line in the same file corresponding to a read and creates a smaller ID string from that line: UMI-Strand-Chromosome-Position```
    Input: [NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC,	0,	2,	76814284,	36,	71M	*,	0,	0,	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA, 6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEEAEEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/, MD:Z:71,	NH:i:1,	HI:i:1,	NM:i:0,	SM:i:36,	XQ:i:40,	X2:i:0,	XO:Z:UU]
    Output: 'CTGTTCAC-0-3-76814284'

# PSEUDOCODE:

1. Make a dictionary which will have the structure : {'UMI-POS-CHROM-STRAND':[whole line from sam file]}

2. Read in a line from the file

3. Store the line as an array with items separated by tabs. Then, for that line:
    a) Check the cigar string:
        - if the cigar string indicates soft clopping (if there is an "S" in the cigar string), edit the position to be the "un-soft-clipped" position: call Correct_Pos()
    b) Create the ID string: 'UMI-Strand-Chromsome-Position' by calling Create_String()
        - get the last 8 characters from item 1 in the read array (UMI), get the strand (item 2), chromosome (item 3), position (item 4)
    c) Check to see if the string already exists as a key in the dictionary
        - If so, continue
        - If not, add dictionary entry for that read ('UMI-POS-CHROM-STRAND-SEQ':[whole line from same file])
    d) Clear the array and move on to the next line
4. At the end of file, write out all the dictionary values






Define the problem
Write examples:
Include a properly formated input sam file
Include a properly formated expected output sam file
Develop your algorithm using pseudocode
Determine high level functions
Description
Function headers
Test examples (for individual functions)
Return statement
For this portion of the assignment, you should design your algorithm for single-end data, with 96 UMIs. UMI information will be in the QNAME, like so: NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT. Discard any UMIs with errors (or error correct, if you're feeling ambitious).