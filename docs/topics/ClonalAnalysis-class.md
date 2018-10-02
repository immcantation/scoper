**ClonalAnalysis-class** - *Output of the analyzeClones function*

Description
--------------------

`ClonalAnalysis` contains output from the [analyzeClones](analyzeClones.md) function.
It includes infromation to interpret clonal assignment performance.




Slots
-------------------



`threshold`
:   cut-off separating the inter (within) and intra (between)
clonal distances.

`inter_intra`
:   data.frame containing all inter and intra clonal distances.

`plot_inter_intra`
:   density plot of inter versus intra clonal distances. The threshold is
shown with a horizental dashed-line.

`neighborhoods`
:   a numeric vector containing scale parameters used in spectral
clustering process.

`plot_neighborhoods`
:   histogram of neighborhoods. The threshold is shown with a vertical
dashed-line.




See also
-------------------

[analyzeClones](analyzeClones.md)



