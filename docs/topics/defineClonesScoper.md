**defineClonesScoper** - *Assigning Ig sequences into clonal groups*

Description
--------------------

The `defineClonesScoper` function provides an unsupervised pipline for assigning Ig sequences into
clonal groups sharing same V gene, J gene, and junction length.


Usage
--------------------
```
defineClonesScoper(db, junction = "JUNCTION", v_call = "V_CALL",
j_call = "J_CALL", first = FALSE, cdr3 = FALSE, mod3 = FALSE,
iter_max = 1000, nstart = 25, nproc = 1, progress = FALSE,
out_name = NULL, out_dir = ".")
```

Arguments
-------------------

db
:   data.frame with Change-O style columns containing sequence data.

junction
:   name of the column containing nucleotide sequences to compare.
Also used to determine sequence length for grouping.

v_call
:   name of the column containing the V-segment allele calls.

j_call
:   name of the column containing the J-segment allele calls.

first
:   if `TRUE` only the first call of the gene assignments
is used. if `FALSE` the union of ambiguous gene
assignments is used to group all sequences with any
overlapping gene calls.

cdr3
:   if `TRUE` remove 3 nts from both ends of `junction`
(converts IMGT junction to CDR3 region). if `TRUE` remove `junction`(s)
with length less than 7 nts.

mod3
:   if `TRUE` remove `junction`(s) with number of nucleotides not modulus of 3.

iter_max
:   the maximum number of iterations allowed for kmean clustering step.

nstart
:   the number of random sets chosen for kmean clustering initialization.

nproc
:   number of cores to distribute the function over.

progress
:   if `TRUE` print a progress bar.

out_name
:   if not `NULL` save cloned data.frame and a summary of cloning
performance. `out_name` string is used as the prefix of the successfully
processed output files.

out_dir
:   specify to change the output directory. The input file
directory is used if this is not specified while `out_name` is specified.




Value
-------------------

Returns a modified `db` data.frame with clone identifiers in the `CLONE` column.
if `out_name` is not `NULL`, it will save the modified `db` and a summary
of cloning performance in the current directory or the specified `out_dir`.


Details
-------------------

An unsupervised pipeline to identify B cell clones from adaptive immune receptor
repertoire sequencing (AIRR-Seq) datasets. This method is based on spectral clustering
of the junction sequences of B cell receptors (BCRs, also referred to as Immunoglobulins,
(Igs)) that share the same V gene, J gene and junction length. It uses an adaptive
threshold that analyzes sequences in a local neighborhood.


Note
-------------------

To assess the performance of clonal assignment process check `analyzeClones`.


References
-------------------


1. coming soon..
 



Examples
-------------------

```R
# clone data using defineClonesScoper function
db <- defineClonesScoper(ExampleDb, junction = "JUNCTION", v_call = "V_CALL",
j_call = "J_CALL", first = TRUE)
```


```
CLONES=  1058
RECORDS=  2000
PASS=  2000
FAIL=  0

```




