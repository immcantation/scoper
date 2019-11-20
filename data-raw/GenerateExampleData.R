# Generate example clones

# Imports
library(alakazam)

#### Generate example database ####

# Load data
ExampleDb <- readChangeoDb("data-raw/ExampleDb.gz")
ExampleDb <- ExampleDb[c("SEQUENCE_ID",
                         "SEQUENCE_IMGT",
                         "GERMLINE_IMGT_D_MASK",
                         "V_CALL",
                         "V_CALL_GENOTYPED",
                         "D_CALL",
                         "J_CALL",
                         "JUNCTION",
                         "JUNCTION_LENGTH",
                         "NP1_LENGTH",
                         "NP2_LENGTH",
                         "SAMPLE",
                         "ISOTYPE",
                         "DUPCOUNT")]

# Save
devtools::use_data(ExampleDb, overwrite=TRUE)

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
