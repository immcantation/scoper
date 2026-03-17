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

#### clone - hierarchicalClones with IUPAC parameter ####

test_that("Test hierarchicalClones with IUPAC parameter", {
    # IUPAC and max_n serve different purposes:
    # - IUPAC: Controls (1) validation (which characters allowed) and 
    #          (2) distance calculation method (fast Hamming vs IUPAC-aware scoring)
    # - max_n: Filters sequences by counting non-ATCG characters using regex "[^ATCG]"
    #          (includes N, ?, and all IUPAC codes - whatever passed validation)
    #
    # Processing order: Validation -> Filtering -> Distance calculation
    #
    # Key use cases tested:
    # 1. IUPAC=FALSE, max_n=0: Strict ATCG-only mode with fast distance calculation
    # 2. IUPAC=TRUE, max_n=0: IUPAC-aware distance (but filters out IUPAC seqs anyway)
    # 3. IUPAC=TRUE, max_n=1+: Main use case - IUPAC codes with proper scoring
    # 4. IUPAC=FALSE, max_n=1+: Allow N/? but reject IUPAC codes (validation blocks them)
    # 5. IUPAC=FALSE/TRUE with standard bases: Backward compatibility (same results)
    
    # Create test data with IUPAC ambiguity codes
    # With threshold=0.15 and length 15, distance > 0.15 requires >=3 mutations (3/15=0.2)
    db_iupac <- data.frame(
        sequence_id = paste0("seq", 1:12),
        v_call = rep("IGHV1-1*01", 12),
        j_call = rep("IGHJ1*01", 12),
        # Use IUPAC ambiguity codes: R(A/G), Y(C/T), W(A/T), S(C/G), M(A/C), K(G/T)
        junction = c(
            # Clone 1: sequences similar to TGTGCAAGCTACTGG (0-1 mutations apart)
            "TGTGCRAGCTACTGG", # R = A or G at pos 5
            "TGTGCRAGCTACTGG", # identical to seq1
            "TGTGCAAGCTACTGG", # standard bases only
            "TGTGCYAGCTACTGG", # Y = C or T at pos 5, 1 mutation from seq3
            "TGTGCYAGCTACTGG", # identical to seq4
            "TGTGCMAGCTACTGG", # M = A or C at pos 5, 1 mutation from seq3
            "TGTGCWAGCTACTGG", # W = A or T at pos 5, 1 mutation from seq3
            
            # Clone 2: sequences similar to ACGTTTGGCCAAACC (3+ mutations from Clone 1)
            "ACGTTTGGCCAAACC", # standard bases, distant from Clone 1
            "ACGTTTGGCCAAACC", # identical to seq8
            "ACGKTTGGCCAAACC", # K = G or T at pos 4, 1 mutation from seq8
            "ACGRTTGGCCAAACC", # R = A or G at pos 4, 1 mutation from seq8
            "ACGTTTSGCCAAACC"  # S = C or G at pos 7, 1 mutation from seq8
        ),
        locus = rep("IGH", 12),
        stringsAsFactors = FALSE
    )

    # Test 1: IUPAC=FALSE should reject IUPAC characters at validation stage
    # Validation (IUPAC=FALSE): Only allows A,T,C,G,N,? - REJECTS IUPAC codes (R,Y,W,S,M,K,etc)
    # Result: Sequences with IUPAC codes fail validation and raise error before any filtering
    # Use case: Strict ATCG-only mode with fast Hamming distance (fastDist_rcpp)
    expect_error(
        hierarchicalClones(db_iupac,
            threshold = 0.15,
            method = "nt", linkage = "single",
            junction = "junction",
            v_call = "v_call", j_call = "j_call",
            IUPAC = FALSE,  # Use fast Hamming distance for ATCG only
            summarize_clones = FALSE
        ),
        "invalid sequence characters"
    )

    # Test 2: IUPAC=TRUE with max_n=0 - validates with IUPAC but filters them out
    # 1. Validation (IUPAC=TRUE): Allows standard bases, N, ?, and IUPAC codes - all sequences pass
    # 2. Filtering (max_n=0): Counts non-ATCG using "[^ATCG]" - removes all sequences with any non-ATCG
    # Result: Only seq3, seq8 and seq9 (standard bases only) remain after filtering
    expect_warning(
        expect_message(
            db_result2 <- hierarchicalClones(db_iupac,
                threshold = 0.15,
                method = "nt", linkage = "single",
                junction = "junction",
                v_call = "v_call", j_call = "j_call",
                IUPAC = TRUE,   # Use IUPAC-aware distance calculation
                max_n = 0,      # Filter out sequences with ANY non-ATCG chars
                summarize_clones = FALSE
            ),
            "Running defineClonesScoper in bulk mode and only keep heavy chains"
        ),
        "Removed 9 sequences with non ATCG characters."
    )

    # Only 3 sequences with standard bases remain
    expect_equal(nrow(db_result2), 3)
    expect_equal(sort(db_result2$sequence_id), c("seq3", "seq8", "seq9"))
    # They are in different clones (distant sequences)
    expect_equal(length(unique(db_result2$clone_id)), 2)

    # Test 3a: Create dataset with N characters to distinguish from IUPAC codes
    db_with_n <- db_iupac
    db_with_n[13,] <- NA
    db_with_n$junction[13] <- "TGTGCNAGCTACTGG"  # Add seq13 with N
    db_with_n$sequence_id <- c(db_with_n$sequence_id[1:12], "seq13")
    db_with_n$v_call <- c(db_with_n$v_call[1:12], "IGHV1-1*01")
    db_with_n$j_call <- c(db_with_n$j_call[1:12], "IGHJ1*01")
    db_with_n$locus <- c(db_with_n$locus[1:12], "IGH")
    
    # Test 3b: IUPAC=TRUE with max_n=1 - allows up to 1 non-ATCG character per sequence
    # This is the main use case: proper IUPAC distance calculation with ambiguity codes
    # 1. Validation (IUPAC=TRUE): Allows standard bases (A,T,C,G) plus N, ?, and IUPAC codes (R,Y,W,S,M,K,etc)
    # 2. Filtering (max_n=1): Counts non-ATCG using "[^ATCG]" - keeps sequences with ≤1 such character
    # 3. Distance: Uses IUPAC-aware scoring (alakazam::pairwiseDist) for proper ambiguity handling
    expect_message(
        db_result3 <- hierarchicalClones(db_with_n,
            threshold = 0.15,
            method = "nt", linkage = "single",
            junction = "junction",
            v_call = "v_call", j_call = "j_call",
            IUPAC = TRUE,   # Use IUPAC-aware distance calculation
            max_n = 1,      # Allow up to 1 non-ATCG character (IUPAC or N)
            summarize_clones = FALSE
        ),
        "Running defineClonesScoper in bulk mode and only keep heavy chains"
    )
    
    # All 13 sequences kept (each has ≤1 non-ATCG character)
    expect_equal(nrow(db_result3), 13)
    
    # Verify that we have 2 different clusters (Clone 1 and Clone 2)
    n_clones <- length(unique(db_result3$clone_id))
    expect_equal(n_clones, 2)
    
    # Identical sequences with IUPAC codes should cluster together
    # Note: rows may be reordered, so reference by sequence_id
    get_clone <- function(seq_id) db_result3$clone_id[db_result3$sequence_id == seq_id]
    
    expect_equal(get_clone("seq1"), get_clone("seq2"))  # seq1=seq2 (both R)
    expect_equal(get_clone("seq4"), get_clone("seq5"))  # seq4=seq5 (both Y)
    expect_equal(get_clone("seq8"), get_clone("seq9"))  # seq8=seq9 (standard)
    
    # seq13 (with N) should cluster with Clone 1 (similar to seq3)
    expect_equal(get_clone("seq3"), get_clone("seq13"))
    
    # Verify two distinct clones based on biological distance
    clone1_id <- get_clone("seq3")  # seq3 in Clone 1
    clone2_id <- get_clone("seq8")  # seq8 in Clone 2
    expect_true(clone1_id != clone2_id)
    expect_true(all(sapply(paste0("seq", 1:7), get_clone) == clone1_id))   # Clone 1: seq1-7
    expect_true(all(sapply(paste0("seq", 8:12), get_clone) == clone2_id))  # Clone 2: seq8-12
    expect_true(get_clone("seq13") == clone1_id)                           # seq13 joins Clone 1

    # Test 4: IUPAC=FALSE with max_n=1+ - allows N/? but rejects other IUPAC codes
    # This demonstrates the two-stage process: validation -> filtering
    # 1. Validation (IUPAC=FALSE): Allows A,T,C,G,N,? but REJECTS IUPAC codes (R,Y,W,S,M,K,etc)
    # 2. Filtering (max_n): Counts non-ATCG using "[^ATCG]" regex on sequences that passed validation
    # Use case: Tolerate low-quality positions (N,?) with fast Hamming distance, but no ambiguity codes
    db_with_n_only <- data.frame(
        sequence_id = paste0("seq", 1:5),
        v_call = rep("IGHV1-1*01", 5),
        j_call = rep("IGHJ1*01", 5),
        junction = c(
            "TGTGCAAGCTACTGG",  # standard bases only
            "TGTGCNAGCTACTGG",  # 1 N (allowed with max_n=1)
            "TGTGC?AGCTACTGG",  # 1 ? (allowed with max_n=1)
            "TGTGCNNGCTACTGG",  # 2 N's (filtered with max_n=1)
            "TGTGCRAGCTACTGG"   # IUPAC R code (should fail validation)
        ),
        locus = rep("IGH", 5),
        stringsAsFactors = FALSE
    )
    
    # Should reject at validation stage due to IUPAC code (R) in seq5
    expect_error(
        hierarchicalClones(db_with_n_only,
            threshold = 0.15,
            method = "nt", linkage = "single",
            junction = "junction",
            v_call = "v_call", j_call = "j_call",
            IUPAC = FALSE,  # Reject IUPAC codes
            max_n = 1,      # Allow up to 1 N or ?
            summarize_clones = FALSE
        ),
        "invalid sequence characters"
    )
    
    # Test with only N/? (no IUPAC codes) - should succeed
    db_n_only_valid <- db_with_n_only[1:4, ]  # Remove seq5 with IUPAC code
    
    expect_warning(
        expect_message(
            db_result5 <- hierarchicalClones(db_n_only_valid,
                threshold = 0.15,
                method = "nt", linkage = "single",
                junction = "junction",
                v_call = "v_call", j_call = "j_call",
                IUPAC = FALSE,  # Use fast Hamming distance
                max_n = 1,      # Allow up to 1 N or ?
                summarize_clones = FALSE
            ),
            "Running defineClonesScoper in bulk mode and only keep heavy chains"
        ),
        "Removed 1 sequences with non ATCG characters."  # seq4 with 2 N's filtered
    )
    
    # Should keep seq1, seq2, seq3 (0-1 non-ATCG characters each)
    expect_equal(nrow(db_result5), 3)
    expect_equal(sort(db_result5$sequence_id), c("seq1", "seq2", "seq3"))
    # All are similar sequences, should be in same clone
    expect_equal(length(unique(db_result5$clone_id)), 1)

    # Test 5: Standard bases work with both IUPAC=TRUE and IUPAC=FALSE
    # This demonstrates backward compatibility and performance consideration:
    # - IUPAC=FALSE: Uses fast Hamming distance (fastDist_rcpp) - faster
    # - IUPAC=TRUE: Uses IUPAC-aware scoring (alakazam::pairwiseDist) - slower but handles ambiguity
    # For standard bases (ATCG), results should be identical
    db_standard <- ExampleDb[1:50, ]

    # Verify db_standard contains only standard bases (not IUPAC ambiguity codes)
    standard_chars <- c("A", "T", "C", "G", "N", "?")
    .is_standard <- function(x) {
        all(unique(strsplit(x, "")[[1]]) %in% standard_chars)
    }
    all_standard <- sapply(db_standard$junction, .is_standard)
    expect_true(all(all_standard),
        info = "ExampleDb contains IUPAC ambiguity codes - update test data or test expectations"
    )

    expect_message(
        db_result_false <- hierarchicalClones(db_standard,
            threshold = 0.15,
            method = "nt", linkage = "single",
            junction = "junction",
            v_call = "v_call", j_call = "j_call",
            IUPAC = FALSE,  # Fast Hamming distance for ATCG only
            summarize_clones = FALSE
        ),
        "Running defineClonesScoper in bulk mode and only keep heavy chains"
    )

    expect_message(
        db_result_true <- hierarchicalClones(db_standard,
            threshold = 0.15,
            method = "nt", linkage = "single",
            junction = "junction",
            v_call = "v_call", j_call = "j_call",
            IUPAC = TRUE,   # IUPAC-aware distance (slower but handles ambiguity)
            summarize_clones = FALSE
        ),
        "Running defineClonesScoper in bulk mode and only keep heavy chains"
    )

    # Results should be identical for standard bases
    # (IUPAC just changes the distance calculation method, not the result)
    expect_equal(nrow(db_result_false), nrow(db_result_true))
    expect_true(all(!is.na(db_result_false$clone_id)))
    expect_true(all(!is.na(db_result_true$clone_id)))
    # Clone assignments should be the same (though IDs may differ)
    expect_equal(length(unique(db_result_false$clone_id)),
                 length(unique(db_result_true$clone_id)))
    
    # Test 6: IUPAC=TRUE with max_n=NULL - no filtering, process all sequences
    # This is the most permissive option for IUPAC data
    # 1. Validation (IUPAC=TRUE): Allows standard bases, N, ?, and all IUPAC codes
    # 2. Filtering (max_n=NULL): No filtering - all sequences kept regardless of character count
    # 3. Distance: Uses IUPAC-aware scoring for proper ambiguity handling
    # Use case: Data with extensive IUPAC codes where you want to process everything
    
    # Create test data with multiple IUPAC codes per sequence
    db_iupac_multi <- data.frame(
        sequence_id = paste0("seq", 1:8),
        v_call = rep("IGHV1-1*01", 8),
        j_call = rep("IGHJ1*01", 8),
        junction = c(
            # Clone 1: sequences with multiple IUPAC codes, similar to each other
            "TGTGCRAGCTRYCTGG",  # R at pos 6, R at pos 12, Y at pos 13 (3 non-ATCG)
            "TGTGCRAGCTRYCTGG",  # identical to seq1
            "TGTGCRAGCTWYNMGG",  # R at pos 6, W at pos 12, Y at pos 13, N at pos 14, M at pos 15 (5 non-ATCG)
            "TGTGCAAGCTACCTGG",  # standard bases only (0 non-ATCG)
            # Clone 2: sequences distant from Clone 1, also with IUPAC codes
            "ACGKTTRGCCNNNYYY",  # K, R, 3 N's, 3 Y's (8 non-ATCG characters!)
            "ACGKTTRGCCNNNYYY",  # identical to seq5
            "ACGGTTAGCCAAACCC",  # standard bases only (0 non-ATCG)
            "ACGKTTSGCCNNNYYY"   # K, S, 3 N's, 3 Y's (8 non-ATCG)
        ),
        locus = rep("IGH", 8),
        stringsAsFactors = FALSE
    )
    
    # With max_n=NULL, all sequences should be kept regardless of IUPAC code count
    expect_message(
        db_result6 <- hierarchicalClones(db_iupac_multi,
            threshold = 0.15,
            method = "nt", linkage = "single",
            junction = "junction",
            v_call = "v_call", j_call = "j_call",
            IUPAC = TRUE,     # Use IUPAC-aware distance calculation
            max_n = NULL,     # No filtering - process all sequences
            summarize_clones = FALSE
        ),
        "Running defineClonesScoper in bulk mode and only keep heavy chains"
    )
    
    # All 8 sequences kept (no filtering applied)
    expect_equal(nrow(db_result6), 8)
    expect_equal(sort(db_result6$sequence_id), paste0("seq", 1:8))
    
    # Should have 2 clones based on biological distance
    expect_equal(length(unique(db_result6$clone_id)), 2)
    
    # Helper function to get clone ID by sequence ID
    get_clone6 <- function(seq_id) db_result6$clone_id[db_result6$sequence_id == seq_id]
    
    # Verify identical sequences cluster together despite many IUPAC codes
    expect_equal(get_clone6("seq1"), get_clone6("seq2"))  # both identical with 3 IUPAC
    expect_equal(get_clone6("seq5"), get_clone6("seq6"))  # both identical with 8 IUPAC
    
    # Verify Clone 1 (seq1-4) and Clone 2 (seq5-8) are distinct
    clone1_id <- get_clone6("seq1")
    clone2_id <- get_clone6("seq5")
    expect_true(clone1_id != clone2_id)
    expect_true(all(sapply(paste0("seq", 1:4), get_clone6) == clone1_id))
    expect_true(all(sapply(paste0("seq", 5:8), get_clone6) == clone2_id))
    
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
        dplyr::mutate(expected_clone_id = dplyr::case_when(
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
        dplyr::mutate(expected_clone_id = dplyr::case_when(
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


