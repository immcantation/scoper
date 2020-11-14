# Load test database
e1 <- new.env()
#load(file.path("tests", "data-tests", "ExampleDb.rda"), envir=e1)
load(file.path("..", "data-tests", "ExampleDb.rda"), envir=e1)
db <- get("ExampleDb", envir=e1)
rm(e1)

# Check for pipelines environment
# pipeline_env <- Sys.getenv("CI") == "true"
# cat("Bitbucket Pipelines:", pipeline_env, "\n")
pipeline_env <- TRUE

#### clone - identicalClones ####

test_that("Test identicalClones", {
    # Truth
    expects <- as.integer(c(20, 21, 23, 26, 27, 28, 30, 44, 50, 100))
    
    # Reproduce example
    db <- identicalClones(ExampleDb, method ="nt", 
                          junction = "junction", v_call = "v_call", 
                          j_call = "j_call", summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        db <- identicalClones(ExampleDb, method ="nt",
                              junction = "junction", v_call = "v_call",
                              j_call = "j_call", summarize_clones = FALSE,
                              nproc=2)
        clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
        expect_identical(clones, expects)
    }
})

#### clone - hierarchicalClones ####

test_that("Test hierarchicalClones", {
    # Truth
    expects <- as.integer(c(7, 8, 8, 8, 8, 9, 10, 11, 12, 683))
    
    # Reproduce example
    db <- hierarchicalClones(ExampleDb, threshold = 0.15,
                             method = "nt", linkage = "single",
                             junction = "junction", 
                             v_call = "v_call", j_call = "j_call", 
                             summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        db <- hierarchicalClones(ExampleDb, threshold = 0.15,
                                 method = "nt", linkage = "single",
                                 junction = "junction",
                                 v_call = "v_call", j_call = "j_call",
                                 summarize_clones = FALSE,
                                 nproc=2)
        clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
        expect_identical(clones, expects)
    }
})

#### clone - spectralClones - novj method ####

test_that("Test spectralClones - novj", {
    # Truth
    expects <- as.integer(c(7, 7, 7, 7, 8, 9, 11, 12, 192, 491))
    
    # Reproduce example
    db <- spectralClones(ExampleDb, method = "novj", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.15,
                         summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    cat(paste(clones))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        db <- spectralClones(ExampleDb, method = "novj",
                             junction = "junction", v_call = "v_call",
                             j_call = "j_call", threshold=0.15,
                             summarize_clones = FALSE,
                             nproc=2)
        clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
        expect_identical(clones, expects)
    }
})

#### clone - spectralClones - vj method ####

test_that("Test spectralClones - vj", {
    # Truth
    expects <- as.integer(c(11, 12, 12, 13, 14, 15, 16, 29, 35, 683))
    
    # Reproduce example
    db <- spectralClones(ExampleDb, method = "vj", 
                         germline = "germline_alignment_d_mask",
                         sequence = "sequence_alignment", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.15,
                         summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    cat(paste(clones))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        db <- spectralClones(ExampleDb, method = "vj",
                             germline = "germline_alignment_d_mask",
                             sequence = "sequence_alignment",
                             junction = "junction", v_call = "v_call",
                             j_call = "j_call", threshold=0.15,
                             summarize_clones = FALSE,
                             nproc=2)
        clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
        expect_identical(clones, expects)
    }
})

