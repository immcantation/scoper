**identicalClones** - *Identical clustering-based method for partitioning Ig sequences into clones.*

Description
--------------------

The `identicalClones` function provides a computational pipline for assigning Ig 
sequences into clonal groups sharing same V gene, J gene, and identical junction.


Usage
--------------------
```
identicalClones(db, method = c("nt", "aa"), junction_col = "junction",
v_call_col = "v_call", j_call_col = "j_call",
clone_col = c("clone_id", "CLONE"), first = FALSE, cdr3 = FALSE,
mod3 = FALSE, max_n = NULL, nproc = 1, verbose = FALSE,
log_verbose = FALSE, out_dir = ".", summarize_clones = FALSE)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

method
:   one of the `"nt"` for nucleotide based clustering or 
`"aa"` for amino acid based clustering.

junction_col
:   character name of the column containing junction sequences.
Also used to determine sequence length for grouping.

v_call_col
:   character name of the column containing the V-segment allele calls.

j_call_col
:   character name of the column containing the J-segment allele calls.

clone_col
:   one of the `"CLONE"` or `"clone_id"` for the output column name 
containing the clone ids.

first
:   specifies how to handle multiple V(D)J assignments for initial grouping. 
If `TRUE` only the first call of the gene assignments is used. 
If `FALSE` the union of ambiguous gene assignments is used to 
group all sequences with any overlapping gene calls.

cdr3
:   if `TRUE` removes 3 nts from both ends of `"junction_col"`
(converts IMGT junction to CDR3 region). if `TRUE` remove 
`junction_col`(s) with length less than 7 nts.

mod3
:   if `TRUE` removes `junction_col`(s) with number of nucleotides not 
modulus of 3.

max_n
:   The maximum number of N's to permit in the junction sequence before excluding the 
record from clonal assignment. Note, under model `"hierarchical"` and method 
`"single"` non-informative positions can create artifactual links between 
unrelated sequences. Use with caution. Default is set to be `"NULL"` for no action.

nproc
:   number of cores to distribute the function over.

verbose
:   if `TRUE` report a summary of each step cloning process;
if `FALSE` process cloning silently.

log_verbose
:   if `TRUE` write verbose logging to a file in `out_dir`.

out_dir
:   specify the output directory to save `log_verbose`. The input 
file directory is used if this is not specified.

summarize_clones
:   if `TRUE` performs a series of analysis to assess the clonal landscape.
See Value for description.




Value
-------------------

For `summarize_clones` = `FALSE`, a modified data.frame with clone identifiers in the `clone_col` column. 
For `summarize_clones` = `TRUE` returns a list containing:

+ `db`:                   modified `db` data.frame with clone identifiers in the `clone_col` column. 
+ `vjl_group_summ`:       data.frame of clones summary, e.g. size, V-gene, J-gene, junction lentgh,
and so on.
+ `inter_intra`:          data.frame containing minimum inter (between) and maximum intra (within) 
clonal distances.
+ `eff_threshold`:        effective cut-off separating the inter (between) and intra (within) clonal 
distances.
+ `plot_inter_intra`:     ggplot histogram of inter (between) versus intra (within) clonal distances. The 
effective threshold is shown with a horizental dashed-line.

If `log_verbose` = `TRUE`, it will write verbose logging to a file in the current directory or 
the specified `out_dir`.


Details
-------------------

`identicalClones` provides a computational platform to explore the B cell clonal 
relationships in high-throughput Adaptive Immune Receptor Repertoire sequencing (AIRR-seq) 
data sets. This function performs clustering among sequences of B cell receptors 
(BCRs, also referred to as Immunoglobulins, (Igs)) that share the same V gene, J gene, and identical junction:



Examples
-------------------

```R
results <- identicalClones(ExampleDb, method ="nt", 
junction_col = "JUNCTION", v_call_col = "V_CALL", 
j_call_col = "J_CALL", summarize_clones = TRUE)
```




