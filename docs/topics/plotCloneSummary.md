**plotCloneSummary** - *Plot ScoperClones object for the summarize_clones=TRUE*

Description
--------------------

`plotCloneSummary` plots the results from `summarize_clones=TRUE` method of 
`spectralClones`, `identicalClones`, and `hierarchicalClones` functions, 
including the minimum inter (between) and maximum intra (within) clonal distances, 
and the calculated efective threshold.


Usage
--------------------
```
plotCloneSummary(
data,
xmin = NULL,
xmax = NULL,
breaks = NULL,
binwidth = NULL,
title = NULL,
size = 0.75,
silent = FALSE,
...
)
```

Arguments
-------------------

data
:   [ScoperClones](ScoperClones-class.md) object output by the [spectralClones](spectralClones.md), 
[identicalClones](identicalClones.md), and [hierarchicalClones](hierarchicalClones.md) functions, 
if argument `summarize_clones` assigned to be `TRUE`.

xmin
:   minimum limit for plotting the x-axis. If `NULL` the limit will 
be set automatically.

xmax
:   maximum limit for plotting the x-axis. If `NULL` the limit will 
be set automatically.

breaks
:   number of breaks to show on the x-axis. If `NULL` the breaks will 
be set automatically.

binwidth
:   binwidth for the histogram. If `NULL` the binwidth 
will be set automatically.

title
:   string defining the plot title.

size
:   numeric value for lines in the plot.

silent
:   if `TRUE` do not draw the plot and just return the ggplot2 
object; if `FALSE` draw the plot.

...
:   additional arguments to pass to ggplot2::theme.




Value
-------------------

A ggplot object defining the plot.



Examples
-------------------

```R
results <- hierarchicalClones(ExampleDb, threshold=0.15,
method="nt", linkage="single",
junction="junction", 
v_call="v_call", j_call="j_call", 
summarize_clones=TRUE)

```

**Error in defineClonesScoper(db = db, threshold = threshold, model = "hierarchical", **: Column(s) 'clone_id' already exist.
 Invalid column names are: 'clone_id', 'VJ_GROUP', 'VJL_GROUP', 'JUNCTION_L', 'CDR3', 'CDR3_L', 'CLONE_temp'.
```R

# Plot clonal summaries 
plot(results, binwidth=0.02)
```

**Error in plot(results, binwidth = 0.02)**: object 'results' not found

See also
-------------------

See [ScoperClones](ScoperClones-class.md) for the the input object definition and 
[spectralClones](spectralClones.md), [identicalClones](identicalClones.md), and [hierarchicalClones](hierarchicalClones.md) 
functions for generating the input object.






