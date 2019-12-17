# CPA2019

C/C++ sequence alignment scoring using Smith-Waterman algorithm.

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
* Gap extension penalty : 1
* BLOSUM matrix : BLOSUM62

The default output will be the top 20 sequence titles with their score.

### Optional command-line parameters

You can customize penalties and BLOSUM matrix on which the algorithm is based.

```
./main database_file query_protein_file_to_base_scores [-o gap_opening_value] [-e gap_extension_value] [-b BLOSUM_matrix_file_path]
```

This line will launch the program with custom gap penalties and BLOSUM matrix. Make sure the
BLOSUM matrix file exists in the given path.

#### All command-line parameters:

- -o gap_opening_value : sets the gap opening value to given argument.
- -e gap_extension_value : sets the gap expansion value to given argument.
- -b BLOSUM_matrix_path : sets the BLOSUM matrix to given path. Needs the matrix file to exist.
- -n N : restricts the database to N sequences.
- -s offset : sets the program to start at given offset in the database.
- -t top : sets the number of top proteins to be displayed.