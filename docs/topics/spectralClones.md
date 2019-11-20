**spectralClones** - *Spectral clustering-based method for partitioning Ig sequences into clones.*

Description
--------------------

The `spectralClones` function provides an unsupervised computational pipline for assigning Ig 
sequences into clonal groups sharing same V gene, J gene, and junction length, based on the 
junction sequence similarity and shared mutations in V and J segments.


Usage
--------------------
```
spectralClones(db, method = c("novj", "vj"),
germline = "germline_alignment", sequence = "sequence_alignment",
junction = "junction", v_call = "v_call", j_call = "j_call",
clone = "clone_id", targeting_model = NULL, len_limit = NULL,
first = FALSE, cdr3 = FALSE, mod3 = FALSE, max_n = NULL,
threshold = NULL, base_sim = 0.95, iter_max = 1000,
nstart = 1000, nproc = 1, verbose = FALSE, log_verbose = FALSE,
out_dir = ".", summarize_clones = FALSE)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

method
:   one of the `"novj"` or `"vj"`. See Details for description.

germline
:   character name of the column containing the germline or reference sequence.

sequence
:   character name of the column containing input sequences.

junction
:   character name of the column containing junction sequences.
Also used to determine sequence length for grouping.

v_call
:   character name of the column containing the V-segment allele calls.

j_call
:   character name of the column containing the J-segment allele calls.

clone
:   the output column name containing the clone ids.

targeting_model
:   [TargetingModel](http://www.rdocumentation.org/packages/shazam/topics/TargetingModel-class) object. Only applicable if `method` = `"vj"`. 
See Details for description.

len_limit
:   [IMGT_V](http://www.rdocumentation.org/packages/shazam/topics/IMGT_SCHEMES) object defining the regions and boundaries of the Ig 
sequences. If NULL, mutations are counted for entire sequence. Only 
applicable if `method` = `"vj"`.

first
:   specifies how to handle multiple V(D)J assignments for initial grouping. 
If `TRUE` only the first call of the gene assignments is used. 
If `FALSE` the union of ambiguous gene assignments is used to 
group all sequences with any overlapping gene calls.

cdr3
:   if `TRUE` removes 3 nts from both ends of `"junction"`
(converts IMGT junction to CDR3 region). if `TRUE` remove 
`junction`(s) with length less than 7 nts.

mod3
:   if `TRUE` removes `junction`(s) with number of nucleotides not 
modulus of 3.

max_n
:   The maximum number of N's to permit in the junction sequence before excluding the 
record from clonal assignment. Note, under model `"hierarchical"` and method 
`"single"` non-informative positions can create artifactual links between 
unrelated sequences. Use with caution. Default is set to be `"NULL"` for no action.

threshold
:   the upper-limit cut-off for clonal grouping.

base_sim
:   required similarity cut-off for sequences in equal distances from each other.

iter_max
:   the maximum number of iterations allowed for kmean clustering step.

nstart
:   the number of random sets chosen for kmean clustering initialization.

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

For `summarize_clones` = `FALSE`, a modified data.frame with clone identifiers in the `clone` column. 
For `summarize_clones` = `TRUE` returns a list containing:

+ `db`:                   modified `db` data.frame with clone identifiers in the `clone` column. 
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

`spectralClones` provides a computational platform to explore the B cell clonal 
relationships in high-throughput Adaptive Immune Receptor Repertoire sequencing (AIRR-seq) 
data sets. Two methods are included to perform clustering among sequences of B cell receptors 
(BCRs, also referred to as Immunoglobulins, (Igs)) that share the same V gene, J gene and junction length: 

+  If `method` = `"novj"`: clonal relationships are inferred using an adaptive threshold that 
indicates the level of similarity among junction sequences in a local neighborhood. 
+  If `method` = `"vj"`: clonal relationships are inferred not only based on the junction region 
homology, but also takes into account the mutation profiles in the V and J segments. Mutation counts are 
determined by comparing the input sequences (in the column specified by `sequence`) to the effective 
germline sequence (IUPAC representation of sequences in the column specified by `germline`). 
+  Not mandatory, but the influence of SHM hot- and cold-spot biases in the clonal inference process will be noted 
if a SHM targeting model is provided through argument `targeting_model` (see [createTargetingModel](http://www.rdocumentation.org/packages/shazam/topics/createTargetingModel) 
for more technical details). 
+  Not mandatory, but the upper-limit cut-off for clonal grouping can be provided to
prevent sequences with disimilarity above the threshold group together.




Examples
-------------------

```R
results <- spectralClones(ExampleDb, method = "vj", 
germline = "germline_alignment_d_mask", 
sequence = "sequence_alignment", 
junction = "junction", v_call = "v_call", 
j_call = "j_call", threshold=0.15, summarize_clones = TRUE)
```




