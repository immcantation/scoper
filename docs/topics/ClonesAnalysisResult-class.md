**ClonesAnalysisResult-class** - *Output of the clonesAnalysis function*

Description
--------------------

`ClonesAnalysisResult` contains output from the [clonesAnalysis](clonesAnalysis.md) function.
It includes infromation to interpret clonal assignment performance.




Slots
-------------------



`threshold`
:   cut-off separating the inter (within) and intra (between)
clonal distances.

`interVsIntra`
:   data.frame containing all inter and intra clonal distances.

`plotInterVsIntra`
:   density plot of inter versus intra clonal distances. The threshold is
shown with a horizental dashed-line.

`neighborhoods`
:   a numeric vector containing scale parameters used in spectral
clustering process.

`plotNeighborhoods`
:   histogram of neighborhoods. The threshold is shown with a vertical
dashed-line.




See also
-------------------

[clonesAnalysis](clonesAnalysis.md)



