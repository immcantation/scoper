# Generate example clones

# Imports
library(alakazam)

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

# Save
use_this::use_data(ExampleDb, overwrite=TRUE)

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
