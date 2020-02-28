# Load test database
e1 <- new.env()
#load(file.path("tests", "data-tests", "ExampleDb.rda"), envir=e1)
load(file.path("..", "data-tests", "ExampleDb.rda"), envir=e1)
db <- get("ExampleDb", envir=e1)
rm(e1)

#ensure older version of sample() used
R_v <- paste(version$major, version$minor,sep=".")
if ( numeric_version(R_v) >= numeric_version("3.6.0") ) {
    RNGkind(sample.kind="Round")   
    set.seed(12345)
}

#### clone - identicalClones ####

test_that("Test identicalClones", {
    ## Reproduce example
    db <- identicalClones(ExampleDb, method ="nt", 
                          junction = "junction", v_call = "v_call", 
                          j_call = "j_call", summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expects <- as.integer(c(20, 21, 23, 26, 27, 28, 30, 44, 50, 100))
    ## Test if the updated function reproduces results
    expect_identical(clones, expects)
})

#### clone - hierarchicalClones ####

test_that("Test hierarchicalClones", {
    ## Reproduce example
    db <- hierarchicalClones(ExampleDb, threshold = 0.15,
                             method = "nt", linkage = "single",
                             junction = "junction", 
                             v_call = "v_call", j_call = "j_call", 
                             summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expects <- as.integer(c(7, 8, 8, 8, 8, 9, 10, 11, 12, 683))
    ## Test if the updated function reproduces results
    expect_identical(clones, expects)
})

#### clone - spectralClones - novj method ####

test_that("Test spectralClones - novj", {
    ## Reproduce example
    db <- spectralClones(ExampleDb, method = "novj", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.15,
                         summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expects <- as.integer(c(7, 7, 7, 7, 8, 9, 10, 12, 192, 491))
    ## Test if the updated function reproduces results
    expect_identical(clones, expects)
})

#### clone - spectralClones - vj method ####

test_that("Test spectralClones - vj", {
    ## Reproduce example
    db <- spectralClones(ExampleDb, method = "vj", 
                         germline = "germline_alignment_d_mask",
                         sequence = "sequence_alignment", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.15,
                         summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expects <- as.integer(c(12, 12, 13, 13, 14, 15, 16, 29, 35, 667))
    ## Test if the updated function reproduces results
    expect_identical(clones, expects)
})


