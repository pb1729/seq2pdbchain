# Make PDB files for unfolded chains given an amino acid sequence

Start from an amino acid sequence and get a corresponding PDB file containing an unfolded protein with that sequence. The program makes an effort to avoid self-collisions and to generate a fairly compact chain. You can view the output and rerun if it is not satisfactory. Currently, the code will still sometimes allow some collisions between side chains and neighbouring amino acids. Running an energy minimization is reccommended!

## Usage:

Pick a destination file. eg. `chains/8P59_unfolded.pdb`. Then run the program:

```
python seq2pdbchain.py > chains/8P59_unfolded.pdb
```

the program will wait for you to paste in an amino acid sequence (use the standard single-letter code) eg:

```
GSHMVPISFVFNRFPRMVRDLAKKMNKEVNFIMRGEDTELDRTFVEEIGEPLLHLLRNAIDHGIEPKEERIAKGKPPIGTLILSARHEGNNVVIEVEDDGRGIDKEKIIRKAIEKGLIDESKAATLSDQEILNFLFVPGFSTKEKVSEVSGRGVGMDVVKNVVESLNGSISIESEKDKGTKVTIRLPLT
```

Now it will run. Progress will be printed to stderr while the actual PDB data will be written to stdout. Be warned, this can take a while. Longer sequences tend to take significantly longer than short ones.

```
current collisions: 65
current collisions: 61
current collisions: 61
current collisions: 50
current collisions: 52
current collisions: 53
current collisions: 53
current collisions: 53
current collisions: 53
current collisions: 49
current collisions: 49
current collisions: 48
current collisions: 47
current collisions: 18
current collisions: 18
current collisions: 16
current collisions: 15
current collisions: 14
current collisions: 14
current collisions: 14
current collisions: 15
current collisions: 14
current collisions: 11
current collisions: 10
current collisions: 10
current collisions: 10
current collisions: 10
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 12
current collisions: 9
current collisions: 9
current collisions: 9
current collisions: 9
current collisions: 9
current collisions: 9
current collisions: 9
current collisions: 9
current collisions: 6
current collisions: 5
current collisions: 5
current collisions: 5
current collisions: 5
current collisions: 5
current collisions: 5
current collisions: 3
current collisions: 3
current collisions: 3
current collisions: 3
current collisions: 3
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 2
current collisions: 1
current collisions: 1
current collisions: 1
current collisions: 0
current radius: 95.527865
current radius: 95.527865
current radius: 89.326070
current radius: 88.637366
current radius: 71.712816
current radius: 65.009812
current radius: 64.518025
current radius: 61.941947
current radius: 60.930471
current radius: 60.556377
current radius: 58.769269
current radius: 58.370336
current radius: 58.333448
current radius: 58.033571
current radius: 55.912477
current radius: 55.494571
current radius: 52.716482
current radius: 52.575870
```




