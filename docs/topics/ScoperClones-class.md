**ScoperClones-class** - *Output of `identicalClones`, `hierarchicalClones`, and `spectralClones`, 
if `summarize_clones=TRUE`*

Description
--------------------

`ScoperClones` contains output from [identicalClones](identicalClones.md), [hierarchicalClones](hierarchicalClones.md), and
[spectralClones](spectralClones.md) functions, if argument `summarize_clones` is `TRUE`.


Usage
--------------------
```
"print"(x)
```
```
"plot"(x, y, ...)
```

Arguments
-------------------

x
:   ScoperClones object

y
:   ignored.

...
:   arguments to pass to [plotCloneSummary](plotCloneSummary.md).




Slots
-------------------



`db`
:   modified input `db` data.frame with clone identifiers in the `clone` 
column.

`vjl_group_summ`
:   data.frame of clones summary, e.g. size, V-gene, J-gene, junction lentgh,
and so on.

`inter_intra`
:   data.frame containing minimum inter (between) and maximum intra (within) 
clonal distances.

`eff_threshold`
:   effective cut-off separating the inter (between) and intra (within) clonal 
distances.




See also
-------------------

[identicalClones](identicalClones.md), [hierarchicalClones](hierarchicalClones.md), and [spectralClones](spectralClones.md)






