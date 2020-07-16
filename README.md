# Compiler for ACDC networks
Example usage: `python compiler.py -i testgraph.txt -n -d 3 -l 15`  
`-i testgraph.txt`: use `testgraph.txt` as input file  
`-n`: run NUPACK locally  
`-d 3`: set NUPACK normalised defect stoppping criterion to 3% (default 5%)
`-l 15`: set central domain length to 15 (default 17)

## Dependencies
- Python 3.7.6
- NumPy 1.18.1
- iGraph 0.7.1
- NUPACK 3.2.2

## Input file structure
Each line in a valid input file contains three symbols from left to right, separated by whitespace:
- The name of some species, say `A`
- An edge symbol, either `->` or `-|`
- The name of another species, say `B`  

Currently only connected networks are supported.
