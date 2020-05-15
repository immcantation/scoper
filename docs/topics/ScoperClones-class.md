**ScoperClones-class** - *S4 class containing clonal assignments and summary data*

Description
--------------------

`ScoperClones` stores output from [identicalClones](identicalClones.md), [hierarchicalClones](hierarchicalClones.md) and
[spectralClones](spectralClones.md) functions.


Usage
--------------------
```
"print"(x)
```
```
"summary"(object)
```
```
"plot"(x, y, ...)
```
```
"as.data.frame"(x)
```

Arguments
-------------------

x
:   ScoperClones object

object
:   ScoperClones object

y
:   ignored.

...
:   arguments to pass to [plotCloneSummary](plotCloneSummary.md).




Slots
-------------------



`db`
:   `data.frame` of repertoire data including with clonal identifiers in 
the column specified during processing.

`vjl_groups`
:   `data.frame` of clonal summary, including sequence count, V gene, 
J gene, junction length, and clone counts.

`inter_intra`
:   `data.frame` containing minimum inter (between) and maximum intra 
(within) clonal distances.

`eff_threshold`
:   effective cut-off separating the inter (between) and intra (within) clonal 
distances.




See also
-------------------

[identicalClones](identicalClones.md), [hierarchicalClones](hierarchicalClones.md) and [spectralClones](spectralClones.md)






