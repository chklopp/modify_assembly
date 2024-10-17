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
save
```

The script can :
- cut : for a given sequence (name) from postion 1 (0 based) to position 2 in the clipboard
- paste : for a give sequence (name) at postion  (0 based) paste clipboard
- reverse_complement : for  a give sequence (name) reverse complement the sequence
- delete : for  a give sequence (name) delete the sequence
- load : load the sequences from the file given using the --input parameter
- save : save the sequences to the file given using the --output parameter
