library(Rcpp)
sourceCpp("src/fastDistAA.cpp")

fastDistAA <- function(seqs) {
  n <- length(seqs)
  v <- fastDistAA_rcpp(seqs)
  structure(v, class="dist", Size=n, Labels=names(seqs), Diag=FALSE, Upper=FALSE)
}

# Naive reference: pure Hamming with X/? semantics matching fastDistAA
# Rules:
#   - same known AA -> match
#   - X vs known AA (either order) -> match (wildcard)
#   - ? vs ? -> match
#   - everything else -> mismatch
naiveDistAA <- function(s1, s2) {
  c1 <- strsplit(toupper(s1), "")[[1]]
  c2 <- strsplit(toupper(s2), "")[[1]]
  known <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")
  mismatches <- 0L
  for (i in seq_along(c1)) {
    a <- c1[i]; b <- c2[i]
    is_match <- (a == b && a %in% known) ||
                (a == "X" && b %in% known) ||
                (b == "X" && a %in% known) ||
                (a == "?" && b == "?")
    if (!is_match) mismatches <- mismatches + 1L
  }
  mismatches
}

# ---- parse args ----
args    <- commandArgs(trailingOnly=TRUE)
verbose <- "-v" %in% args
k_flag  <- which(args == "-k")
l_flag  <- which(args == "-l")
K       <- if (length(k_flag) && k_flag < length(args)) as.integer(args[k_flag + 1L]) else 500L
L       <- if (length(l_flag) && l_flag < length(args)) as.integer(args[l_flag + 1L]) else 20L

cat(sprintf("Testing fastDistAA: k=%d sequences, l=%d length\n", K, L))

set.seed(42)
AAS  <- c("A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y","X","?")
seqs <- replicate(K, paste(sample(AAS, L, replace=TRUE), collapse=""))

# ---- pairwise naive reference ----
pairs <- combn(K, 2)  # 2 x choose(K,2) matrix
naive <- integer(ncol(pairs))
for (p in seq_len(ncol(pairs))) {
  naive[p] <- naiveDistAA(seqs[pairs[1,p]], seqs[pairs[2,p]])
}

# ---- fastDistAA output (lower-triangle, column-major) ----
fast <- as.integer(fastDistAA(seqs))

# combn gives column-major lower-triangle in the same order R's dist uses
errors <- 0L
for (p in seq_len(ncol(pairs))) {
  i <- pairs[1,p]; j <- pairs[2,p]
  if (verbose)
    cat(sprintf("%s  %s  naive=%d  fast=%d\n",
                seqs[i], seqs[j], naive[p], fast[p]))
  if (fast[p] != naive[p]) {
    cat(sprintf("FAIL pair (%d,%d): '%s' vs '%s'  fast=%d  naive=%d\n",
                i, j, seqs[i], seqs[j], fast[p], naive[p]))
    errors <- errors + 1L
  }
}

if (errors == 0L) {
  cat(sprintf("All %d pairs passed.\n", ncol(pairs)))
} else {
  cat(sprintf("%d / %d pairs FAILED.\n", errors, ncol(pairs)))
  quit(status=1)
}
