# CRISPR_distance
An implementation of the Kupczok and Bollback (2013) method.

Kupczok, A., and Bollback, J.P. 2013. Probabilistic models for CRISPR spacer content evolution. BMC Evolutionary Biology 13(1): 54. doi:10.1186/1471-2148-13-54.

Reference paper:

Perez, M., Angers B., Young C.R., and Juniper S.K. (2021) Shining light on a deep-sea bacterial symbiont population structure with CRISPR. Microbial Genomics 7(8): 000625. doi:10.1099/mgen.0.000625

## DEPENDENCIES:
python 3.6
python libraries:
pandas
numpy
biopython (Bio.Phylo)
scipy
mpmath
itertools
math
sys
time
argparse
multiprocessing

## INSTALLATION:
```
# Clone this repository
$ git clone https://github.com/maepz/CRISPR_distance.git

# Go into the repository
$ cd CRISPR_distance
```

## RUNNING THE PROGRAM:

```
python Run_Parallel_CRISPR_Distance.py Exemples/test_arrays.txt Exemples/Test_minimal
```
