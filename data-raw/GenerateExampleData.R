# Generate example clones

# Imports
library(alakazam)
library(dplyr)

#### Generate example database ####

# Load data
ExampleDb <- readChangeoDb("data-raw/ExampleDb.gz")
ExampleDb <- ExampleDb[c("sequence_id",
                         "sequence_alignment",
                         "germline_alignment_d_mask",
                         "v_call",
                         "v_call_genotyped",
                         "d_call",
                         "j_call",
                         "junction",
                         "junction_length",
                         "np1_length",
                         "np2_length",
                         "sample_id",
                         "c_call",
                         "duplicate_count")]
ExampleDb$locus <- getLocus(ExampleDb$v_call)

c_trans <- c(IGHM="IgM", IGHD="IgD", IGHA="IgA", IGHG="IgG")
ExampleDb <- ExampleDb %>%
    mutate(c_call=translateStrings(c_call, c_trans),
           germline_alignment=germline_alignment_d_mask)


# Save
usethis::use_data(ExampleDb, overwrite=TRUE)

#### Generate example clones ####
# ClonedExampleDb <- defineClonesScoper(db,
#                                       junction = "JUNCTION",
#                                       v_call = "V_CALL",
#                                       j_call = "J_CALL",
#                                       first = TRUE,
#                                       cdr3 = FALSE,
#                                       mod3 = FALSE,
#                                       nproc = 1,
#                                       progress = T)
# # Save
# devtools::use_data(ClonedExampleDb, overwrite=TRUE)
