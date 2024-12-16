## Simple Script for parsing and analyzing TRAVIS Voronoi analsis
It keeps the molecular types and indicies from the voro.txt that list all the voronoi cells in contact with the voronoi cell of the reference atom. The shells of atoms from the same reference molecule are united and listed in a DataFrame. 
For each molecular type (the ones detected by travis) 2 csv file are created: <name>_full.csv and <name>_size.csv; first list all shells with their TRAVIS index, the second list each shell composition and their occurrence in a one-shot encoding where rows are reference molecule and columns are observed molecules in the coordination shell

###Usage
```
import Analyzer from shell
obj = Analyzer("voro.txt")
```
