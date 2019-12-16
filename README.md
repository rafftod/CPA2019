# CPA2019

C++ sequence alignment scoring using Smith-Waterman algorithm.

## Compiling

After cloning the repository, run the following command :

```
make
```

The program is now compiled and ready to run.

## Running the program

```
./main database_file query_protein_file_to_base_scores
```

This will run the program with default options:
* Gap open penalty : 11
* Gap expansion penalty : 1
* BLOSUM matrix : BLOSUM62

The output will be the top 20 sequence titles with their score.

### Optional command-line parameters

You can customize penalties and BLOSUM matrix on which the algorithm is based.

```
./main database_file query_protein_file_to_base_scores [-o gap_opening_value] [-e gap_expansion_value] [-b BLOSUM_matrix_file_path]
```

This line will launch the program with custom gap penalties and BLOSUM matrix. Make sure the
BLOSUM matrix file exists in the given path.