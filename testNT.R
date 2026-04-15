library(Rcpp)
sourceCpp("src/fastDist.cpp")

fastDist <- function(seqs) {
  n <- length(seqs)
  v <- fastDist_rcpp(seqs)
  structure(v, class="dist", Size=n, Labels=names(seqs), Diag=FALSE, Upper=FALSE)
}

# Naive reference implementation
trueDist <- function(s1, s2) {
  c1 <- strsplit(s1, "")[[1]]
  c2 <- strsplit(s2, "")[[1]]
  known <- c("A", "C", "G", "T")
  mismatches <- 0L
  for (i in seq_along(c1)) {
    a <- c1[i]; b <- c2[i]
    match <- (a == b && a %in% known) ||
             (a == "N" && b %in% known) ||
             (b == "N" && a %in% known) ||
             (a == "?" && b == "?")
    if (!match) mismatches <- mismatches + 1L
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

cat(sprintf("Testing fastDist: k=%d sequences, l=%d length\n", K, L))

set.seed(42)
BASES <- c("A", "C", "G", "T", "N", "?")
seqs  <- replicate(K, paste(sample(BASES, L, replace=TRUE), collapse=""))

# ---- pairwise naive reference ----
pairs <- combn(K, 2)
naive <- integer(ncol(pairs))
for (p in seq_len(ncol(pairs))) {
  naive[p] <- trueDist(seqs[pairs[1,p]], seqs[pairs[2,p]])
}

# ---- fastDist output ----
fast <- as.integer(fastDist(seqs))

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
