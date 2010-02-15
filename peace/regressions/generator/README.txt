input_generator.py: 
Generates a file of regression sequences based on the generation of
blocks (randomly generated sequences of a fixed size) and a template
on how to put the blocks together.

Requires biopython to run.

Parameters:
* Block length (presumably window size)
* Input file name (template)
* Output file name (fasta file)

Template format:
Tuples defining block sequences (one per line)
* Each tuple is a sequence of integers (comma separated), specifying 
  the blocks that will be concatenated together
* Block numbers are indexed from 0
* A negative entry specifies the use of the reverse complement of a
block
* Comments can be inserted by starting the line with #
