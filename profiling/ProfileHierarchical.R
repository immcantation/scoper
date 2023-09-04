# Imports
library(scoper)
library(profvis)

#### Load example data ####
# setwd("profiling")
data <- readr::read_tsv("1000.tsv", col_types = readr::cols())
data$LOCUS <- alakazam::getLocus(data$V_CALL)

#### Calculate expected mutations ####
profvis({
    hierarchicalClones(data, threshold=0.12, junction = "JUNCTION",
                       v_call = "V_CALL", j_call = "J_CALL", locus = "LOCUS",
                       linkage = "single", normalize = "len", nproc = 1)
})

