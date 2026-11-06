# modify_assembly
simple script to cut, paste, reverse complement, delete contigs in an assembly

When you compare different genome assemblies you can spot differences which are due to scaffolding errors rather the read chromosomal rearrangement. This script has been developped to ease large scaffold or contig modification. It takes a fasta input file and command file as inputs and outputs a modified fasta output file.

```python manipulate_sequence.py --input seq1a.fasta --output seq1a.out.fasta --commands commands_example.txt```

The command file contains a set of commands which are played by the script in the top down order 

```
load seq1a.fasta
cut seq1 2 20
paste seq2 0
reverse_complement seq1
cut seq1 10 40
paste seq2 50
clean seq2
filter 100
save
```

The script can :
- cut : for a given sequence (name) from postion 1 (0 based) to position 2 in the clipboard
- paste : for a give sequence (name) at postion  (0 based) paste clipboard
- reverse_complement : for  a give sequence (name) reverse complement the sequence
- delete : for  a give sequence (name) delete the sequence
- load : load the sequences from the file given using the --input parameter
- split : for a given sequence (name) and position (0-based), split the sequence in two; the first keeps the original name and the second gets -2 suffix added to the name
- move : for a given sequence (name) and 1-based position, move the sequence to that position in the file (1 = first position)
- rename : for a given sequence (old name) and new name, rename the sequence while preserving its position in the file
- filter : for a given lower limit (size), remove all sequences smaller than this size from the file
- join : for two given sequences (first name, second name) and the resulting name, join the second sequence to the first one; the two source sequences are removed and the resulting sequence takes the position of the first one
- save : save the sequences to the file given using the --output parameter
