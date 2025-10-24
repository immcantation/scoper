# Load test database
e1 <- new.env()

load(file.path("..", "data-tests", "ExampleDb.rda"), envir=e1)
load(file.path("..", "data-tests", "Example10x.rda"), envir=e1)
db <- get("ExampleDb", envir=e1)
db_sc <- get("Example10x", envir=e1)
rm(e1)

#ensure older version of sample() used
R_v <- paste(version$major, version$minor,sep=".")
if ( numeric_version(R_v) >= numeric_version("3.6.0") ) {
    expect_warning(RNGkind(sample.kind="Round"))
}

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
    expects <- c(6, 7, 7, 7, 7, 8, 10, 11, 12, 679)
    
    # Reproduce example
    set.seed(12345)
    db <- spectralClones(ExampleDb, method = "novj", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.20,
                         summarize_clones = FALSE)
    clones  <- as.numeric(tail(sort(table(db$clone_id)), 10))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        set.seed(12345)
        db <- spectralClones(ExampleDb, method = "novj",
                             junction = "junction", v_call = "v_call",
                             j_call = "j_call", threshold=0.20,
                             summarize_clones = FALSE,
                             nproc=2)
        clones  <- as.numeric(tail(sort(table(db$clone_id)), 10))
        expect_identical(clones, expects)
    }
})

#### clone - spectralClones - vj method ####

test_that("Test spectralClones - vj", {
    # Truth
    expects <- as.integer(c(11, 12, 12, 13, 14, 15, 16, 29, 35, 683))
    
    # Reproduce example
    set.seed(12345)
    db <- spectralClones(ExampleDb, method = "vj", 
                         germline = "germline_alignment_d_mask",
                         sequence = "sequence_alignment", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", threshold=0.15,
                         summarize_clones = FALSE)
    clones <- as.integer(as.vector(tail(sort(table(db$clone_id)), 10)))
    expect_identical(clones, expects)
    
    # Test parallel
    if (!pipeline_env) {
        set.seed(12345)
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

#### Single cell 

test_that("Test assigning clones works for heavy-only sc data", {
    # Issue https://bitbucket.org/kleinstein/scoper/issues/23
    db_sc$chain <- "light"
    db_sc$chain[grepl("IGH",db_sc[['v_call']])] <- "heavy"
    db_sc_heavy <- db_sc %>%
        filter(chain == "heavy")
    expect_warning(cloned <- identicalClones(db_sc_heavy, method="aa",
                               cell_id = "cell_id",
                               locus = "locus", nproc=1),
                   "Single cell mode requested, but")
})


#### Single cell 

test_that("Test prepare_db", {
    
    # Note:
    # seq7 (cell4) has v_call IGHV3*01,IGHV1*01 and will be assigned G2 when first=F. 
    #     The order is not what we would usually expect, because it is not sorted
    #     alpahbetically. It is not sorted on purpose just to test the first= T/F parameter.
    db <- data.frame(
        sequence_id=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8"),
        cell_id=c("cell1","cell1","cell2","cell2","cell3","cell3","cell4","cell4"),
        v_call=c("IGHV1*01","IGLV1*01","IGHV1*01","IGLV1*01","IGHV1*01","IGLV2*01", "IGHV3*01,IGHV1*01","IGLV1*01"),
        d_call=c("IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA),
        j_call=c("IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01*01","IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01"),
        junction=c("TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTTTTC","TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTC")
    )
    db$chain <- "light"
    db$chain[grepl("IGH",db[['v_call']])] <- "heavy"
    db$locus <- alakazam::getLocus(db$v_call)
    db$junction_len <- stringi::stri_length(db[['junction']])

    # Case: heavy T, first F
    # prepare_db uses groupGenes to find vj gene groups, but
    # uses junc_len = NULL
    # It adds the junction length groups outside groupGenes
    # requires both cell_id and locus
    # scoper:::prepare_db
    db_only_heavy_T_locus_first_F <- scoper:::prepare_db(
        db,
        only_heavy = T,
        cell_id = "cell_id",
        locus="locus", 
        first=F)$db
    expect_true(all(db_only_heavy_T_locus_first_F$vj_group =="G1"))
    
    # Case: heavy T, first T
    # scoper:::prepare_db

    db_only_heavy_T_locus_first_T <- scoper:::prepare_db(
        db,
        only_heavy = T,
        cell_id = "cell_id",
        locus = "locus",
        first = T
    )$db

    # The warning
    # "no non-missing arguments to max; returning -Inf"
    # comes from an igraph function inside groupGenes.
    # It is not expected, but can be ignored. It could depend on the
    # version of igraph installed (seen with 2.2.0).
    # Example:
    # r$> igraph::graph_from_adjacency_matrix(adjmatrix = mtx_adj,
    #             mode = "undirected", diag = FALSE )
    # IGRAPH 7ecf7be UN-- 2 0 --
    # + attr: name (v/c)
    # + edges from 7ecf7be (vertex names):
    # Warning message:
    # In max(el[, 3]) : no non-missing arguments to max; returning -Inf
    #
    # It seems to happen when mtx_adj is and identity matrix (only 1s in the 
    # diagonal and no edges) and diag=FALSE (no self-links, the diagonal will
    # be zeroed out and max, used somewhere in the function, returns and error).
    # r$> mtx_adj
    # 2 x 2 sparse Matrix of class "dgCMatrix"
    #             IGHV1@IGHJ1 IGHV3@IGHJ1
    # IGHV1@IGHJ1           1           .
    # IGHV3@IGHJ1           .           1
    #
    # r$> max()
    # [1] -Inf
    # Warning message:
    # In max() : no non-missing arguments to max; returning -Inf

    expect_equal(db_only_heavy_T_locus_first_T$vj_group, c("G1","G1","G1","G1","G1","G1","G2","G2"))


    # Case: heavy F, first T
    # Expectation: throw an error as only_heavy=F is not supported any more.

    expect_warning(
        db_only_heavy_F_locus_first_T <- scoper:::prepare_db(
            db,
            only_heavy = F,
            cell_id = "cell_id",
            locus = "locus",
            first = T
        )$db,
        "The only_heavy = FALSE parameter is deprecated. Will run as if only_heavy = TRUE"
    )

    expect_equal(db_only_heavy_F_locus_first_T$vj_group, c("G1","G1","G1","G1","G1","G1","G2","G2"))
    
    # Case: heavy F, first F
    # these are also all G1.....
    expect_warning(
        db_only_heavy_F_locus_first_F <- scoper:::prepare_db(db,
                                                       only_heavy = F, 
                                                       cell_id="cell_id",
                                                       locus="locus",
                                                       first=F)$db,
        "The only_heavy = FALSE parameter is deprecated. Will run as if only_heavy = TRUE"
    )

    expect_equal(db_only_heavy_F_locus_first_F$vj_group, c("G1","G1","G1","G1","G1","G1","G1","G1"))
    
    # Results heavy T/F with same first setting should be the same, as heavy is ignored
    expect_equal(db_only_heavy_T_locus_first_T$vj_group, db_only_heavy_F_locus_first_T$vj_group)
    expect_equal(db_only_heavy_T_locus_first_F$vj_group, db_only_heavy_F_locus_first_F$vj_group)
    
})

test_that("Test hierarchicalClones only_heavy and first", {
    # https://testthat.r-lib.org/articles/third-edition.html
    # https://testthat.r-lib.org/articles/third-edition.html#warnings
    # https://testthat.r-lib.org/articles/third-edition.html#messages
    local_edition(3)

    # Note:
    # seq7 (cell4) has v_call IGHV3*01,IGHV1*01 and will be assigned G2 when first=F. 
    #     The order is not what we would usually expect, because it is not sorted
    #     alpahbetically. It is not sorted on purpose just to test the first= T/F parameter.
    db <- data.frame(
        sequence_id=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8"),
        cell_id=c("cell1","cell1","cell2","cell2","cell3","cell3","cell4","cell4"),
        v_call=c("IGHV1*01","IGLV1*01","IGHV1*01","IGLV1*01","IGHV1*01","IGLV2*01", "IGHV3*01,IGHV1*01","IGLV1*01"),
        d_call=c("IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA),
        j_call=c("IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01*01","IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01"),
        junction=c("TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTTTTC","TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTC")
    )
    db$chain <- "light"
    db$chain[grepl("IGH",db[['v_call']])] <- "heavy"
    db$locus <- alakazam::getLocus(db$v_call)
    db$junction_len <- stringi::stri_length(db[['junction']])

    ## Test hierachicalClones
    # CGJ 1/29/25 updated case to remove split_light testing
    # split_light = FALSE is set to avoid the warning message
    ## case only_heavy   first
    ## 1    T            F    
    ## 2    T            T    
    ## 3    F            F   
    ## 4    F            T    

    # Case 1
    expect_message(
        clones <- hierarchicalClones(
            db,
            threshold = 0,
            cell_id = "cell_id",
            locus = "locus",
            only_heavy = TRUE,
            split_light = FALSE,
            summarize_clones = TRUE,
            first = F, # default is first=F
            nproc = 1
        ),
        "Running defineClonesScoper in single cell mode",
        fixed = TRUE
    )
    expect_true(all(clones@db[['clone_id']] == "1"))

    # Case 2
    expect_message(
        clones <- hierarchicalClones(
            db,
            threshold = 0,
            cell_id = "cell_id",
            locus = "locus",
            only_heavy = TRUE,
            split_light = FALSE,
            summarize_clones = TRUE,
            first = T, # default is first=F
            nproc = 1
        ),
        "Running defineClonesScoper in single cell mode",
        fixed = TRUE
    )
    expect_equal(clones@db[["clone_id"]], c("1", "1", "1", "1", "1", "1", "2", "2"))

    # Case 3
    expect_message(
        expect_warning(
            expect_warning(
                clones <- hierarchicalClones(
                    db,
                    threshold = 0,
                    cell_id = "cell_id",
                    locus = "locus",
                    only_heavy = FALSE,
                    split_light = TRUE,
                    summarize_clones = TRUE,
                    first = F, # default is first=F
                    nproc = 1
                ),
                "split_light = TRUE is deprecated. Please use split_light = FALSE. After clonal identification, light chain groups can be found with dowser::resolveLightChains"
            ),
            "only_heavy = FALSE is deprecated. Running as if only_heavy = TRUE"
        ),
        "Running defineClonesScoper in single cell mode",
    )

    expect_true(all(clones@db[["clone_id"]] == "1"))
    
    # Case 4
    expect_warning(
        expect_warning(
            expect_message(
                clones <- hierarchicalClones(
                    db,
                    threshold = 0,
                    cell_id = "cell_id",
                    locus = "locus",
                    only_heavy = FALSE,
                    split_light = TRUE,
                    summarize_clones = TRUE,
                    first = T, # default is first=F
                    nproc = 1
                ),
                "Running defineClonesScoper in single cell mode",
                fixed = TRUE
            ),
            "only_heavy = FALSE is deprecated. Running as if only_heavy = TRUE"
        ),
        "split_light = TRUE is deprecated. Please use split_light = FALSE. After clonal identification, light chain groups can be found with dowser::resolveLightChains"
    )

    expect_equal(clones@db[['clone_id']],c("1","1","1","1","1","1","2","2"))
    
    #TODO: check in dowser if comment is addressed.
    ## What happens with the second groupGenes (inside the light chain split) when
    # first=false, but the "linker" ambiguous call was left out of the same cluster id
    # because of the distance threshold? This groupGenes could be splitting again by vj calls...
    # seq2 is the linker, and the juction is one nt different. Everything else is the same.
    db <- data.frame(
        sequence_id=c("seq1","seq2","seq3","seq4","seq5","seq6"),
        cell_id=c("1","2","3","1","2","3"),
        v_call=c("IGHV1", "IGHV1,IGHV2","IGHV2","IGLV1","IGLV1","IGLV2"),
        j_call=c("IGHJ1","IGHJ1","IGHJ1","IGHJ1","IGHJ1","IGHJ1"),
        junction=c("TCGAAATTC","TCGAACTTC","TCGAAATTC","TCGAAATTC","TCGAAATTC","TCGAAATTC")
    )
    db$locus <- alakazam::getLocus(db$v_call)
    db

    # If threshold = 0, seq1 (cell 1) and seq3 (cell 3) are in the same clone, seq2 is separate.
    expect_message(
        clones_split_F_th0 <- hierarchicalClones(
            db,
            threshold = 0,
            cell_id = "cell_id",
            locus = "locus",
            only_heavy = TRUE,
            split_light = FALSE,
            summarize_clones = TRUE,
            first = F, # default is first=F
            nproc = 1
        )
    )
    
    # expect sequences from cells 1 and 3 in clone 1, and
    # sequences from cell 2 in clone 2
    clones_split_F_th0@db %>%
        mutate(expected_clone_id = dplyr::case_when(
            cell_id %in% c(1, 3) ~ "1",
            cell_id %in% c(2) ~ "2",
            TRUE ~ NA_character_
        )) %>%
        {
            expect_equal(.$clone_id, .$expected_clone_id)
        }

    # If threshold = 0.12, all sequences are in the same clone
    expect_message(
        clones_split_F_th012 <- hierarchicalClones(
            db,
            threshold = 0.12,
            cell_id = "cell_id",
            locus = "locus",
            only_heavy = TRUE,
            split_light = FALSE,
            summarize_clones = TRUE,
            first = F, # default is first=F
            nproc = 1
        )
    )
    expect_true(all(clones_split_F_th012@db[['clone_id']] == "1"))

    expect_message(
        expect_warning(
            clones_split_T_th0 <- hierarchicalClones(
                db,
                threshold = 0,
                cell_id = "cell_id",
                locus = "locus",
                only_heavy = TRUE,
                split_light = TRUE,
                summarize_clones = TRUE,
                first = F, # default is first=F
                nproc = 1
            ),
            "split_light = TRUE is deprecated. Please use split_light = FALSE. After clonal identification, light chain groups can be found with dowser::resolveLightChains"
        ),
        "Running defineClonesScoper in single cell mode"
    )
    # expecting same results because split_ligth not supported anymore.
    expect_equal(clones_split_F_th0@db[['clone_id']], clones_split_T_th0@db[['clone_id']])
    
    ## TODO: multiple light chains?
    db2 <- bind_rows(
        db,
        data.frame( # adding a second light chain for cell 3
            sequence_id = c("seq6"),
            cell_id = c("3"),
            v_call = c("IGLV1,IGLV2"),
            j_call = c("IGHJ1"),
            junction = c("TCGAAATTC"),
            locus="IGL"
        )
    )

   # If threshold = 0, seq1 (cell 1) and seq3 (cell 3) are in the same clone, seq2 is separate.
    expect_message(
        clones_split_F_th0_2l <- hierarchicalClones(
            db2,
            threshold = 0,
            cell_id = "cell_id",
            locus = "locus",
            only_heavy = TRUE,
            split_light = FALSE,
            summarize_clones = TRUE,
            first = F, # default is first=F
            nproc = 1
        )
    )
    
    # expect sequences from cells 1 and 3 in clone 1, and
    # sequences from cell 2 in clone 2, regardless of
    # light chains
    clones_split_F_th0_2l@db %>%
        mutate(expected_clone_id = dplyr::case_when(
            cell_id %in% c(1, 3) ~ "1",
            cell_id %in% c(2) ~ "2",
            TRUE ~ NA_character_
        )) %>%
        {
            expect_equal(.$clone_id, .$expected_clone_id)
        }
})

## Add test for clones by light chain
# CGJ 1/29/25 This is no longer needed as we do not split by light chains so I 
# changed the test to best for a warning when split_lights = TRUE
test_that("Testing split_light warnings for all cloning mehtods", {
  # Load test db with light chain ambiguity
  db <- data.frame(
    sequence_id=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8"),
    cell_id=c("cell1","cell1","cell2","cell2","cell3","cell3","cell4","cell4"),
    v_call=c("IGHV1*01","IGLV1*01","IGHV1*01","IGLV1*01","IGHV1*01","IGLV2*01", "IGHV3*01,IGHV1*01","IGLV1*01"),
    d_call=c("IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA,"IGHD1*01",NA),
    j_call=c("IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01*01","IGHJ1*01","IGLJ1*01","IGHJ1*01","IGLJ1*01"),
    junction=c("TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTTTTC","TCGAAATTC","TCGTTTTTC","TCGAAATTC","TCGTTTTTC")
  )
  db$chain <- "light"
  db$chain[grepl("IGH",db[['v_call']])] <- "heavy"
  db$locus <- alakazam::getLocus(db$v_call)
  db$junction_len <- stringi::stri_length(db[['junction']])  
  # Run hierarchicalClones without and with splitting light chain
   db_nsplit <- hierarchicalClones(
      db,
      threshold=0.1,
      cell_id='cell_id',
      locus='locus',
      only_heavy=TRUE,
      split_light=FALSE,
      summarize_clones=TRUE,
      first=F, # default is first=F
      nproc=1)
   expect_warning(db_split <- hierarchicalClones(
     db,
     threshold=0.1,
     cell_id='cell_id',
     locus='locus',
     only_heavy=TRUE,
     split_light=TRUE,
     summarize_clones=TRUE,
     first=F, # default is first=F
     nproc=1))
   # make sure they are the same 
   expect_equal(db_nsplit@db[['clone_id']], db_split@db[['clone_id']])
   
   # Run spectralClones with and without splitting light chain
   set.seed(12345)
   db_nsplit <- spectralClones(db, method = "novj", 
                        junction = "junction", v_call = "v_call", 
                        j_call = "j_call", threshold=0.20, cell_id = "cell_id",
                        summarize_clones = FALSE, split_light = FALSE)
   expect_warning(db_split <- spectralClones(db, method = "novj", 
                               junction = "junction", v_call = "v_call", 
                               j_call = "j_call", threshold=0.20, cell_id = "cell_id",
                               summarize_clones = FALSE, split_light = TRUE))
   expect_equal(db_nsplit[['clone_id']], db_split[['clone_id']])
   
   # Run identicalClones with and withouth splitting light chain 
   db_nsplit <- identicalClones(db, method ="nt", 
                         junction = "junction", v_call = "v_call", 
                         j_call = "j_call", summarize_clones = FALSE, 
                         cell_id = "cell_id", split_light = FALSE)
   expect_warning(db_split <- identicalClones(db, method ="nt", 
                                junction = "junction", v_call = "v_call", 
                                j_call = "j_call", summarize_clones = FALSE, 
                                cell_id = "cell_id", split_light = TRUE))
   expect_equal(db_nsplit[['clone_id']], db_split[['clone_id']])
})
