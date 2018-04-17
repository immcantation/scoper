library("methods")
library("doParallel")

setClass("InterVsInteraResult",
         slots = c(interVsIntera="list",
                   threshold="numeric",
                   meanInter="numeric",
                   sdInter="numeric",
                   meanIntera="numeric",
                   sdIntera="numeric"))

setClass("clonalAnalysisResult",
         slots = c(interVsIntera="list",
                   threshold="numeric",
                   meanInter="numeric",
                   sdInter="numeric",
                   meanIntera="numeric",
                   sdIntera="numeric",
                   plotInterVsIntera="list",
                   neighborhoods="numeric",
                   plotNeighborhoods="list"))

#################################
# Define universal plot settings
#################################
getBaseTheme <- function() {
    base_theme <- theme_bw() + 
        theme(text=element_text(size=5),
              plot.margin=grid::unit(c(0.5,0.8,0.5,0.5), units="line")) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              panel.border = element_rect(size = 2.0, linetype = "solid", colour = "black", fill = NA)) +
        theme(strip.background=element_blank(),
              strip.text=element_text(size=16, face="bold")) +
        theme(legend.position='bottom',
              legend.spacing=grid::unit(2, "points"),
              legend.text=element_text(size=12),
              legend.title=element_text(size=15),
              legend.key.height=grid::unit(10, "points"), 
              legend.key.width=grid::unit(20, "points"))
    return(base_theme)
}

#################################
# sub-functions 
#################################
filterFunction <- function(groupBy, filterBy) {
    paste(paste(groupBy, paste("'",filterBy,"'", sep=""), sep="=="), collapse=" & ")
}

#################################
# epsilon calculator by "null"
#################################
epsilonNullD <- function(d, vec) {
    return(sd(vec))
}

#################################
# epsilon calculator by "assign"
#################################
assign <- function(d, vec) {
    vec <- vec[vec<=d]
    std <- ifelse (length(vec) == 1, 0, sd(vec)) 
    return(std)
}

#################################
# epsilon calculator by "infer"
#################################
infer <- function(d, vec) {
    vec <- sort(vec)
    diffVec <- diff(vec)
    if (length(unique(diffVec[diffVec > 0])) == 1) {
        d <- ceiling(mean(vec))
    } else {
        id <- which.max(diffVec)   
        d <- vec[id]
    }
    vec <- vec[vec<=d]
    std <- ifelse (length(vec) == 1, 0, sd(vec)) 
    return(std)
}

#################################
# kernel matrix calculator
#################################
krnlMtxGenerator <- function(mtx, neighborhood=c("assign", "infer"), d=NULL) {
    # Radial basis function kernel: In the Gaussian Kernel if two points are 
    # close then K_ij≈1 and when two points are far apart then Kij≈0
    fun <- match.arg(neighborhood)
    n <- dim(mtx)[1]
    # calculate epsilons
    epsilon <- rep(0, length=n)
    for (i in 1:n) {
        epsilon[i] <- apply(as.matrix(mtx[i,]), 2 , FUN=fun, d=d)  
        #as.matrix produces a single-column matrix, therefore we use 2 in apply.  
    } 
    # calculate kernel matrix
    krnl_mtx <- matrix(data=1, nrow=n, ncol=n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            krnl_mtx[i,j] <- exp(-mtx[i,j]^2/(2*epsilon[i]*epsilon[j]))
            krnl_mtx[j,i] <- krnl_mtx[i,j]
        }
    }
    # if mtx[i,j] and epsilon == 0
    krnl_mtx[is.nan(krnl_mtx)] <- 1   
    # Disconnect those edges with distance larger than d
    if (neighborhood == "assign") krnl_mtx[mtx > d] <- 0
    return(krnl_mtx)
}

#################################
 # affinity matrix calculator 
#################################
makeAffinity <- function(mtx, n.neighboors=2) {
    n <- nrow(mtx)
    if (n.neighboors >= n) {  # fully connected
        aff_mtx <- mtx
    } else {
        aff_mtx <- matrix(data=0, ncol=n, nrow=n)
        for(i in 1:n) { 
            # for each line only connect to those points with larger similarity 
            best.similarities <- sort(mtx[i,], decreasing=TRUE)[1:n.neighboors]
            for (s in best.similarities) {
                j <- which(mtx[i,] == s)
                aff_mtx[i,j] <- mtx[i,j]
                aff_mtx[j,i] <- mtx[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
            }
        }
    }
    return(aff_mtx)
}

#################################
# laplacian matrix calculator
#################################
laplacian_mtx <- function(entry) {
    # Calculate unnormalised Laplacian matrix and its eigenfunctions
    D <- diag(apply(entry, 1, sum))
    L <- D - entry
    return(L)
}

#################################
# order of magnitude calculator
#################################
log10_ceiling <- function(x) {
    ceiling(log10(x))-1
}

#################################
# lowest integer
#################################
floor_dec <- function(x, level=1) {
    round(x - 5*10^(-level-1), level)
}

#################################
# assigne v and j jenes 
#################################
Parse_VJ <- function(db, 
                     vCallColumn="V_CALL", 
                     jCallColumn="J_CALL", 
                     first=TRUE) {
    if (first) {
        db$V <- alakazam::getGene(db[[vCallColumn]])
        db$J <- alakazam::getGene(db[[jCallColumn]])
    } else {
        db$V1 <- alakazam::getGene(db[[vCallColumn]], first=FALSE)
        db$J1 <- alakazam::getGene(db[[jCallColumn]], first=FALSE)
        db$V <- db$V1
        db$J <- db$J1
        # Reassign V genes to most general group of genes
        for(ambig in unique(db$V1[grepl(',', db$V1)])) {
            for(g in strsplit(ambig, split=',')[[1]]) {
                db$V[grepl(g, db$V1)] = ambig
            }
        }
        # Reassign J genes to most general group of genes
        for(ambig in unique(db$J1[grepl(',',db$J1)])) {
            for(g in strsplit(ambig, split=',')[[1]]) {
                db$J[grepl(g, db$J1)] = ambig
            }
        }
    }
    return(db)
}

#################################
# main cloning functions
#################################
defineCloneSpectral <- function(db, 
                                junctionColumn = "JUNCTION", 
                                vCallColumn = "V_CALL",
                                jCallColumn = "J_CALL",
                                fields = NULL, 
                                first = TRUE, 
                                cdr3 = FALSE, 
                                mod3 = TRUE,
                                iter.max = 1000, 
                                nstart = 25,
                                nproc = 1, 
                                progress = FALSE, 
                                verbose = FALSE,
                                pathToLog = NULL) {

    # set Seed for reproducibility
    set.seed(12345)
    
    # number of enteries
    rec_count <- nrow(db)
    
    # number of fails
    fail_count <- 0
    
    # Initial checks
    if (all(c(progress, verbose))) { 
        stop("Both progressBar and verbose cannot be TRUE") 
        }
    if (nproc > 1 & all(verbose == TRUE | !is.null(pathToLog))) {  
        stop("if nproc > 1, only progress bar can be seen (progress=TRUE, verbose=FALSE, pathToLog=NULL).") 
    }
    if (nproc > 1 & verbose == TRUE) {
        stop("if nproc > 1, only progress bar can be seen (progress=TRUE, verbose=FALSE, pathToLog=NULL).") 
    }
    if (!is.data.frame(db)) { 
        stop("Must submit a data frame") 
    }
    if ("CLONE" %in% colnames(db)) { 
        stop("Column 'CLONE' already exist.") 
    }
    
    # neighborhood <- match.arg(neighborhood)
    neighborhood <- "infer"
    similarity <- NULL
    # if (is.null(neighborhood)) { 
    #     stop(" 'neighborhood' must be specified.") 
    # }
    # if (neighborhood == "assign" & is.null(similarity)) stop("similarity needs to be assigned.")
    
    if (!is.null(pathToLog)) {
        if (!dir.exists(pathToLog)) { stop("directory", pathToLog, "does not exist") }
        log_path <- file.path(pathToLog, "LOG.txt")
        cat("",
            file=log_path, append=F)
    }
    
    # Check for valid columns
    columns <- c(junctionColumn, vCallColumn, jCallColumn, fields)
    columns <- columns[!is.null(columns)]
    check <- shazam:::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }

    # Check for invalid characters
    valid_seq <- sapply(db[[junctionColumn]], shazam:::allValidChars, colnames(getDNAMatrix(gap=0)))
    not_valid_seq <- which(!valid_seq)
    if (length(not_valid_seq)>0) {
        stop("Invalid sequence characters in the ", junctionColumn, " column. ",
                length(not_valid_seq)," sequence(s) found.", "\n Valid characters are: '",  colnames(getDNAMatrix(gap=0)), "'")
    } 
    
    # Convert sequence columns to uppercase
    db <- shazam:::toupperColumns(db, junctionColumn)  

    # add junction temp column
    db$JUNC_temp <- db[[junctionColumn]]
    
    # add junction length column
    db$L <- stringi::stri_length(db$JUNC_temp)

    # check for mod3 
    if (mod3) {
        db <- filter(db, L%%3 == 0)
    }

    # check for cdr3 
    if (cdr3) {
        # filter junctions with length >6
        db <- filter(db, L > 6)
        # add cdr3 temp column
        db$JUNC_temp <- substr(db$JUNC_temp, 4, db$L-3)
        # update cdr3 length column
        db$L <- stringi::stri_length(db$JUNC_temp)
    }
    
    # Parse V and J columns to get gene
    db <- Parse_VJ(db, 
                   vCallColumn=vCallColumn, 
                   jCallColumn=jCallColumn, 
                   first=first)

    # groups to use
    groupBy <- c("V", "J", "L")
    if (!is.null(fields)) {
        groupBy <- append(groupBy, fields)
    }
    
    # unique groups
    uniqueGroups <- data.frame(unique(db[, groupBy]))
    colnames(uniqueGroups) <- groupBy
    rownames(uniqueGroups) <- NULL
    n_groups <- nrow(uniqueGroups)
    groupBy_length <- length(groupBy)
    
    # Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    # check the progressbar
    if (progress) {
        cat("CLONING> ", "\n")
        pb <- shazam:::progressBar(n_groups)
    }
    # cloning
    db_cloned <- foreach(i=1:n_groups, 
                         .combine="rbind", 
                         .export=c("filterFunction", "spectralClustering", "krnlMtxGenerator", 
                                   "makeAffinity", "kMeanClustering", "laplacian_mtx", "floor_dec", neighborhood),
                         .packages = c("alakazam","plyr","dplyr"), 
                         .errorhandling='stop') %dopar% {
                             ft <- uniqueGroups[i, ]
                             filterBy <- c()
                             for (j in 1:groupBy_length) filterBy <- c(filterBy, as.character(ft[[groupBy[j]]]))
                             db_group <- db %>% 
                                 dplyr::filter_(filterFunction(groupBy, filterBy))
                             db_group$ID <- db_group %>%
                                 dplyr::group_by(JUNC_temp) %>%
                                 dplyr::group_indices()
                             if (length(unique(db_group$JUNC_temp)) == 1) {
                                 CLONE <- data.frame(CLONE=as.vector(paste(i, rep(1, times=nrow(db_group)), sep="_")))
                             } else {
                                 idCluster <- spectralClustering(entrySeq=db_group$JUNC_temp, 
                                                                 id=db_group$ID,
                                                                 similarity=similarity, 
                                                                 neighborhood=neighborhood,
                                                                 iter.max=iter.max, 
                                                                 nstart=nstart)
                                 CLONE <- data.frame(CLONE=as.vector(paste(i, idCluster, sep="_")))
                             }
                             # Update progress
                             if (progress) { pb$tick() }
                             # print verboose
                             if (verbose) {
                                 cat("=================================", "\n", sep="")
                                 cat("group ", i, "/", n_groups, "\n", sep="")
                                 cat("vGene=", uniqueGroups[i,1], "\n", sep="")
                                 cat("jGene=", uniqueGroups[i,2], "\n", sep="")
                                 cat("Junction length=", uniqueGroups[i,3], "\n", sep="")
                                 cat("Number of sequence=", nrow(db_group), "\n", sep="")
                                 cat("Number of unique sequence=", length(unique(db_group$JUNC_temp)), "\n", sep="")
                                 cat("Number of clone=", length(unique(CLONE$CLONE)), "\n", sep="")
                             }
                             # print into log
                             if (!is.null(pathToLog)) {
                                 cat("=================================",
                                     file=log_path, append=T, sep="\n")
                                 cat(paste0("group ", i, "/", n_groups),
                                     file=log_path, append=T, sep="\n")
                                 cat(paste("vGene=", uniqueGroups[i,1]),
                                     file=log_path, append=T, sep="\n")
                                 cat(paste("jGene=", uniqueGroups[i,2]),
                                     file=log_path, append=T, sep="\n")
                                 cat(paste("Junction length=", uniqueGroups[i,3]),
                                     file=log_path, append=T, sep="\n")
                                 cat(paste("Number of sequence=", nrow(db_group)),
                                     file=log_path, append=T, sep="\n")
                                 cat(paste("Number of unique sequence=", length(unique(db_group$JUNC_temp))),
                                     file=log_path, append=T, sep="\n")
                                 cat(paste("Number of clone=", length(unique(CLONE$CLONE))),
                                     file=log_path, append=T, sep="\n")
                             }
                             # return result from each proc
                             return(bind_cols(db_group, CLONE))
                         }
    
    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }

    db_cloned$CLONE_temp <- db_cloned %>%
        dplyr::group_by(CLONE) %>%
        dplyr::group_indices()
    db_cloned$CLONE <- db_cloned$CLONE_temp
    db_cloned <- db_cloned[order(db_cloned$CLONE), ]
    db_cloned$CLONE <- as.character(db_cloned$CLONE)
    
    clone_count <- length(unique(db_cloned$CLONE))
    pass_count <- nrow(db_cloned)
    fail_count <- rec_count - pass_count
    
    cat(paste("CLONES= ", clone_count), "\n", sep="")
    cat(paste("RECORDS= ", rec_count), "\n", sep="")
    cat(paste("PASS= ", pass_count), "\n", sep="")
    cat(paste("FAIL= ", fail_count), "\n", sep="")
    return(db_cloned[, !(names(db_cloned) %in% c("V", "J", "L", "JUNC_temp", "NUMBRE_OF_N", "NUMBRE_OF_DOT", "V1", "J1", "ID","CLONE_temp"))])
}

#################################
# spectral clustering algorithm
#################################
spectralClustering <- function(entrySeq, 
                               id,
                               similarity = NULL, 
                               neighborhood = c("assign", "infer"),
                               iter.max = 25, 
                               nstart = 25) {
    # constants
    nSeq <- length(entrySeq)
    uniqueSeq <- unique(entrySeq)
    nUniqueSeq <- length(uniqueSeq)
    l <- unique(nchar(uniqueSeq))
    uniqueId <-  unique(id)
    # scale parameteres
    neighborhood <- match.arg(neighborhood)
    if (neighborhood == "assign") {
        if (is.null(similarity)) stop("similarity needs to be assigned.")
        d <- round((1-similarity)*l) 
    } else if (neighborhood == "infer") {
        d <- NULL
    }
    # calculate distance matrix
    dist_mtx <- pairwiseDist(seq=uniqueSeq, dist_mat=getDNAMatrix(gap=0))
    # calculate kernel matrix
    krnl_mtx <- krnlMtxGenerator(mtx=dist_mtx, neighborhood=neighborhood, d=d)
    # calculate affinity matrix. n.neighboors could be ceiling(sqrt(nUniqueSeq))
    aff_mtx <- makeAffinity(mtx=round(krnl_mtx, 5), n.neighboors=nUniqueSeq)
    # clustering
    unique_idCluster <- kMeanClustering(entry=aff_mtx, iter.max=iter.max, nstart=nstart)
    # back to reality
    idCluster <- rep(NA, nSeq)
    for (i in 1:nUniqueSeq) {
        j <- which(id == uniqueId[i])
        idCluster[j] <- unique_idCluster[i]
    }
    # return results
    return(idCluster)
}

#################################
# kmean clustering algorithm
#################################
kMeanClustering <- function (entry, 
                             iter.max = 25, 
                             nstart = 25) {
    n <- nrow(entry)
    if (all(entry[!diag(nrow(entry))] == 0)) {
        # affinity matrix is diagonal. Each sequence belongs to a singlton clone.
        return(c(1:n))
    } else {
        # every real symmetric matrix is Hermitian, and therefore all its eigenvalues are real.
        L <- laplacian_mtx(entry)
        evL <- eigen(L, symmetric=TRUE)
        eigenVecs <- round(evL$vectors, 5)
        eigenVals <- floor_dec(evL$values, 3)
        egvRev <- rev(eigenVals)
        # nonZero.ID <- which(egvRev > 0)
        # k <- min(nonZero.ID)-1
        k <- which.max(diff(egvRev))   # slope: forward differencing
        kentry <- eigenVecs[ ,(n-k+1):n]
        clust <- kmeans(kentry, centers=k, iter.max=iter.max, nstart=nstart)
        return(clust$cluster)
    }
}

#################################
# plot eigenvalues function
#################################
plotEigenValuesSpectrum <- function(db, 
                                    vGene, jGene, junctionLength,
                                    first = TRUE,
                                    cdr3 = FALSE) {
    
    # Parse V and J columns to get gene
    db <- Parse_VJ(db,
                   vCallColumn="V_CALL",
                   jCallColumn="J_CALL",
                   first=first)
    
    # filter group
    db_group <- db %>% 
        dplyr::filter(V == vGene, J == jGene, L == junctionLength)
    entrySeq <- db_group$JUNCTION
    # constants
    nSeq <- length(entrySeq)
    uniqueSeq <- unique(entrySeq)
    nUniqueSeq <- length(uniqueSeq)
    if (nUniqueSeq == 1) {
        message(paste("number of unique sequences = 1"))
        return(NULL)
    }
    l <- unique(nchar(uniqueSeq))
    # cdr3
    if (cdr3) uniqueSeq <- substr(uniqueSeq, 4, nchar(uniqueSeq)-3)
    # scale parameteres
    # neighborhood <- match.arg(neighborhood)
    neighborhood <- "infer"
    similarity <- NULL
    if (neighborhood == "assign") {
        if (is.null(similarity)) stop("similarity needs to be assigned.")
        d <- round((1-similarity)*l) 
    } else if (neighborhood == "infer") {
        d <- NULL
    }
    # calculate distance matrix
    dist_mtx <- pairwiseDist(seq=uniqueSeq, dist_mat=getDNAMatrix(gap=0))
    # calculate kernel matrix
    krnl_mtx <- krnlMtxGenerator(mtx=dist_mtx, neighborhood=neighborhood, d=d)
    # calculate affinity matrix. n.neighboors could be ceiling(sqrt(nUniqueSeq))
    aff_mtx <- makeAffinity(mtx=round(krnl_mtx, 5), n.neighboors=nUniqueSeq)
    n <- nrow(aff_mtx)
    if (all(aff_mtx[!diag(nrow(aff_mtx))] == 0)) {
        k <- n
        message(paste("Affinity matrix is diagonal. Unique sequences belong to separate clones. nClone=", k))
        return(list(k, "no eigenvalue", db_group))
    } else {
        L <- laplacian_mtx(aff_mtx)
        evL <- eigen(L, symmetric=TRUE)
        eigenVecs <- round(evL$vectors, 5)
        eigenVals <- floor_dec(evL$values, 3)
        egvRev <- rev(eigenVals)
        # nonZero.ID <- which(egvRev > 0)
        # k <- min(nonZero.ID)-1
        k <- which.max(diff(egvRev))   # slope: forward differencing
        df <- data.frame(id=c(1:n), egv=egvRev)
        p <- ggplot(df, aes(x=id, y=egv)) +
            theme_bw() +
            ggtitle(paste0("nClone=", k)) +
            xlab("Number of Cluster") +
            ylab("EigenValue") +
            geom_point(shape=1, size=1.0) +
            geom_vline(xintercept=k, color="red", linetype=2)
        return(list(p, rev(eigenVals), db_group))
    }
}

#################################
# inter-clone-distance vs intera-clone-distance
#################################
calculateInterVsIntera <- function(db, 
                                   junctionColumn = "JUNCTION", 
                                   vCallColumn = "V_CALL",
                                   jCallColumn = "J_CALL",
                                   first = TRUE,
                                   cdr3 = FALSE,
                                   nproc = 1, 
                                   progress = FALSE, 
                                   verbose = FALSE) {
    # groups to use
    groupBy <- c("V", "J", "L")
    
    # find unique groups of sequences with same vgene, jgene, and junction length
    uniqueGroups <- db %>%
        dplyr::group_by(V, J, L) %>%
        dplyr::summarise(CLONE=paste(unique(CLONE), collapse=","))
    n_groups <- nrow(uniqueGroups)
    
    # Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ 
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    # check the progressbar
    if (progress) {
        cat("INTER AND INTERA DISTANCES ANALYSIS> ", "\n")
        pb <- shazam:::progressBar(n_groups)
    }
    
    # open dataframes
    vec_ff <- foreach(k=1:n_groups,
                      .combine="c",
                      .packages = c("alakazam","plyr","dplyr"), 
                      .errorhandling='stop') %dopar% {
                          clones <- strsplit(uniqueGroups$CLONE[k], split=",")[[1]]
                          l <- uniqueGroups$L[k]
                          n_clones <- length(clones)
                          seqs <- db$JUNC_temp[db$CLONE %in% clones]
                          names(seqs) <- db$CLONE[db$CLONE %in% clones]
                          if (verbose) cat("group=", k, "/", n_groups, ", number of clones=", n_clones, ", number of sequences=", length(seqs), "\n", sep="")
                          seqs_db <- data_frame(value = seqs, name = names(seqs)) %>%
                              group_by(name, value) %>% # alternatively: group_by(name) if name value pair is always unique
                              slice(1) %>%
                              ungroup()
                          seqs <- seqs_db$value 
                          names(seqs) <- seqs_db$name
                          dist_mtx <- pairwiseDist(seqs, dist_mat=getDNAMatrix(gap=0))
                          # prealoocate a dataframe
                          nrow_f <- n_clones + n_clones*(n_clones-1)/2
                          vec_f <- rep(NA, nrow_f)
                          # minimum and maximum distance in each clone
                          n <- 0
                          if (n_clones == 1) {
                              if (!all(dist_mtx == 0)) {
                                  n <- n+1
                                  vec_f[n] <- max(dist_mtx)/l 
                                  names(vec_f)[n] <- paste(clones[1], "inter", sep="_")
                              }
                          } else {
                              for (i in 1:(n_clones-1)) {
                                  xx <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[i]]
                                  if (!all(xx == 0)) {
                                      n <- n+1
                                      vec_f[n] <- max(xx)/l 
                                      names(vec_f)[n] <- paste(clones[i], "inter", sep="_")
                                  }
                                  for (j in (i+1):n_clones) {
                                      xy <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[j]]
                                      xy <- xy[xy>0]
                                      if (length(xy) != 0) {
                                          n <- n+1
                                          vec_f[n] <- min(xy)/l
                                          names(vec_f)[n] <- paste(clones[i], clones[j], "intera", sep="_")
                                      }
                                  }
                              }
                              yy <- dist_mtx[rownames(dist_mtx) == clones[j], colnames(dist_mtx) == clones[j]]
                              if (!all(yy == 0)) {
                                  n <- n+1
                                  vec_f[n] <- max(yy)/l 
                                  names(vec_f)[n] <- paste(clones[j], "inter", sep="_")
                              }   
                          }
                          # Update progress
                          if (progress) { pb$tick() }
                          # remove all na's
                          vec_f <- vec_f[!is.na(vec_f)]
                          # return result from each proc
                          if (length(vec_f) == 0) return(NULL)
                          return(vec_f)
                      }
    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    # convert to a data.frame
    db_dff <- data.frame(keyName=names(vec_ff), VALUE=vec_ff, row.names=NULL)
    db_dff$LABEL <- "inter"
    db_dff$LABEL[grepl("intera", db_dff$keyName)] <- "intera"
    db_dff <- cbind(stringr::str_split_fixed(db_dff$keyName, "_", n=3), db_dff)
    db_dff$keyName <- NULL
    db_dff$`3` <- NULL
    db_dff <- db_dff %>%
        dplyr::rename(CLONE_X=`1`,
                      CLONE_Y=`2`)
    db_dff$CLONE_Y[grepl("inter", db_dff$CLONE_Y)] <- NA
    interVsIntera <- list()
    interVsIntera[[length(interVsIntera)+1]] <- db_dff
    
    # find threshold
    db.summ <- db_dff %>%
        dplyr::group_by(LABEL) %>%
        dplyr::summarise(MEAN=mean(VALUE, na.rm = TRUE),
                         SD=sd(VALUE, na.rm = TRUE))
    func1.1 <- db.summ$MEAN[db.summ$LABEL == "inter"]
    func1.2 <- db.summ$SD[db.summ$LABEL == "inter"]
    func2.1 <- db.summ$MEAN[db.summ$LABEL == "intera"]
    func2.2 <- db.summ$SD[db.summ$LABEL == "intera"]
    minInt <- 0
    maxInt <- 1
    intxn <- uniroot(shazam:::intersectPoint, interval = c(minInt, maxInt), tol=1e-8, extendInt="yes",
                     first_curve = "norm", second_curve = "norm", 
                     func1.0=1, func1.1=func1.1, func1.2=func1.2, 
                     func2.0=1, func2.1=func2.1, func2.2=func2.2)
    threshold <- round(intxn$root, 2)
    meanInter <- round(func1.1, 2)
    sdInter <- round(func1.2, 2)
    meanIntera <- round(func2.1, 2)
    sdIntera <- round(func2.2, 2)
    
    InterVsInteraResult<-new("InterVsInteraResult",
                             interVsIntera=interVsIntera,
                             threshold=threshold,
                             meanInter=meanInter,
                             sdInter=sdInter,
                             meanIntera=meanIntera,
                             sdIntera=sdIntera)
    return(InterVsInteraResult)    
}

#################################
# plot inter-clone-distance vs intera-clone-distance 
#################################
plotInterVsIntera <- function(db, threshold, meanInter, sdInter, meanIntera, sdIntera, 
                              hline_size=0.5, txt_size=3.0) {
    
    # Check for valid columns
    columns <- c("VALUE", "LABEL")
    columns <- columns[!is.null(columns)]
    check <- shazam:::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    db <- select(db, c("VALUE", "LABEL"))
    
    pdat <- db %>%
        dplyr::group_by(LABEL) %>%
        do(data.frame(loc = density(.$VALUE, na.rm = TRUE)$x,
                      dens = density(.$VALUE, na.rm = TRUE)$y
        )
        )
    normalize <- function(x) {
        return ((x-min(x)) / (max(x)-min(x)))
    }
    pdat <- pdat %>%
        dplyr::group_by(LABEL) %>%
        dplyr::mutate(dens_norm=normalize(dens))
    # Flip and offset densities for the groups
    pdat$dens_norm <- ifelse(pdat$LABEL == 'inter', pdat$dens_norm*(-1), pdat$dens_norm)
    # fill color
    fill_manual <- c("inter"="grey30",
                     "intera"="grey60")
    # labels
    mg <- "Mean \u00b1 SD"
    Encoding(mg) <- "UTF-8"
    pm <-"\u00b1"
    Encoding(pm)<-"UTF-8"
    mg1_txt <- data.frame(X=-0.9, Y=threshold-0.025, LAB=paste(mg))
    db1_txt <- data.frame(X=-0.9, Y=threshold-0.05, LAB=paste(meanInter, pm, sdInter))
    mg2_txt <- data.frame(X=0.9, Y=threshold+0.05, LAB=paste(mg))
    db2_txt <- data.frame(X=0.9,  Y=threshold+0.025, LAB=paste(meanIntera, pm, sdIntera))
    
    # plot
    p <- ggplot(pdat) + 
        getBaseTheme() + 
        theme(axis.title.x = element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y=element_text(size=12),
              axis.text.y=element_text(size=11)) +
        ylab("Normalized hamming distance") +
        scale_fill_manual(name="",
                          values=fill_manual, 
                          labels = c("inter"="maximum-distance\nwithin clones  ", 
                                     "intera"="minimum-distance\nbetween clones  ")) + 
        geom_polygon(aes(x=dens_norm, y=loc, fill=LABEL)) +
        geom_hline(yintercept=threshold, color="grey30", linetype=2, size=hline_size) +
        geom_text(data=mg1_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE) +
        geom_text(data=db1_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE) +
        geom_text(data=mg2_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE) +
        geom_text(data=db2_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE)
    return(p)    
}

#################################
# neighborhoods
#################################
calculateNeighborhoods <- function(db, 
                                   junctionColumn = "JUNCTION", 
                                   vCallColumn = "V_CALL",
                                   jCallColumn = "J_CALL",
                                   similarity = NULL, 
                                   neighborhood = c("assign", "infer"),
                                   first = TRUE,
                                   cdr3 = FALSE,
                                   progress = FALSE,
                                   verbose = FALSE) {
    # groups to use
    groupBy <- c("V", "J", "L")
    # find unique groups of sequences with same vgene, jgene, and junction length
    uniqueGroups <- data.frame(unique(db[, groupBy]))
    colnames(uniqueGroups) <- groupBy
    rownames(uniqueGroups) <- NULL
    n_groups <- nrow(uniqueGroups)
    
    # scale parameteres
    neighborhood <- match.arg(neighborhood)
    if (neighborhood == "assign") {
        if (is.null(similarity)) stop("similarity needs to be assigned.")
        d <- round((1-similarity)*l) 
    } else if (neighborhood == "infer") {
        d <- NULL
    }
    
    # check the progressbar
    if (progress) {
        cat("\n")
        cat("NEIGHBORHOOD ANALYSIS> ", "\n")
        pb <- shazam:::progressBar(n_groups)
    }
    
    # sigma anlysis
    eps <- rep(NA, nrow(db)) 
    idx <- 1
    for (i in 1:n_groups) {
        if (verbose) print(paste(n_groups, i, sep=" "))
        db_group <- db %>% 
            dplyr::filter(V == uniqueGroups$V[i],
                          J == uniqueGroups$J[i],
                          L == uniqueGroups$L[i])
        if (nrow(db_group) == 1) {
            # Update progress
            if (progress) { pb$tick() }
            next
        }
        l <- unique(db_group$L)
        entrySeq <- db_group$JUNC_temp
        mtx <- pairwiseDist(seq=entrySeq, dist_mat=getDNAMatrix(gap=0))
        n <- dim(mtx)[1]
        epsilon <- rep(0, length=n)
        for (j in 1:n) {
            epsilon[j] <- apply(as.matrix(mtx[j,]), 2 , FUN=neighborhood, d=d)  
        }
        eps[idx:(idx+n-1)] <- epsilon/l
        idx <- idx+n
        # Update progress
        if (progress) { pb$tick() }
    }
    # remove all na's
    eps <- eps[!is.na(eps)]
    return(eps)
 }

#################################
# plot neighborhoods
#################################
plotNeighborhoods <- function(sigmas, threshold = NULL, vline_size = 0.75, binwidth=0.01, center=0.005) {
    sigmas_df <- data.frame(sigmas=sigmas)
    p <- ggplot(sigmas_df) +
        getBaseTheme() +
        theme(axis.text=element_text(size=11),
              axis.title=element_text(size=12)) +
        xlab("Normalized hamming distance") +
        ylab("Density") +
        geom_histogram(aes(x=sigmas, y=..density..),
                       binwidth=binwidth, center=center, alpha=1.0, fill="black", color="white", size=0.25)
    if (!is.null(threshold)) p <- p + geom_vline(xintercept=threshold, linetype=2, size=vline_size)
    return(p)
}

#################################
# clonal analysis
#################################
clonalAnalysis <- function(db,
                           junctionColumn = "JUNCTION", 
                           vCallColumn = "V_CALL",
                           jCallColumn = "J_CALL",
                           first = TRUE,
                           cdr3 = FALSE,
                           nproc = 1, 
                           progress = FALSE, 
                           verbose = FALSE,
                           hline_size = 0.75, 
                           vline_size = 0.75,
                           txt_size = 3.0,
                           binwidth=0.01, 
                           center=0.005) {
    
    # Initial checks
    if (all(c(progress, verbose))) { 
        stop("Both progressBar and verbose cannot be TRUE") 
    }
    if (nproc > 1 & verbose == TRUE) {
        stop("if nproc > 1, only progress bar can be seen (progress=TRUE, verbose=FALSE).") 
    }
    
    # Check for valid columns
    clone <- "CLONE"
    columns <- c(junctionColumn, vCallColumn, jCallColumn, clone)
    columns <- columns[!is.null(columns)]
    check <- shazam:::checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    
    # add junction temp column
    db$JUNC_temp <- db[[junctionColumn]]
    
    # add junction length column
    db$L <- stringi::stri_length(db$JUNC_temp)
    
    # check for cdr3 
    if (cdr3) {
        # add cdr3 temp column
        db$JUNC_temp <- substr(db$JUNC_temp, 4, db$L-3)
        # update cdr3 length column
        db$L <- stringi::stri_length(db$JUNC_temp)
    }
    
    # Parse V and J columns to get gene
    db <- Parse_VJ(db,
                   vCallColumn=vCallColumn,
                   jCallColumn=jCallColumn,
                   first=first)
    
    # calculate inter and intera distances
    results <- calculateInterVsIntera(db,
                                      junctionColumn = junctionColumn, 
                                      vCallColumn = vCallColumn,
                                      jCallColumn = jCallColumn,
                                      first = first,
                                      cdr3 = cdr3,
                                      nproc = nproc, 
                                      progress = progress, 
                                      verbose = verbose)
    
    # revoke the results
    df <- results@interVsIntera[[1]]
    threshold <- results@threshold
    meanInter <- results@meanInter
    sdInter <- results@sdInter
    meanIntera <- results@meanIntera
    sdIntera <- results@sdIntera
    interVsIntera <- list()
    interVsIntera[[length(interVsIntera)+1]] <- df
    
    # plot inter and intera distances
    p1 <- list()
    p <- plotInterVsIntera(df,
                           threshold=threshold,
                           meanInter=meanInter,
                           sdInter=sdInter,
                           meanIntera=meanIntera,
                           sdIntera=sdIntera ,
                           hline_size=hline_size, 
                           txt_size=txt_size)
    p1[[length(p1)+1]] <- p
    
    # calculate neighborhoods
    neighborhood <- "infer"
    similarity <- NULL
    neighborhoods <- calculateNeighborhoods(db, 
                                            junctionColumn = junctionColumn, 
                                            vCallColumn = vCallColumn,
                                            jCallColumn = jCallColumn,
                                            similarity = similarity, 
                                            neighborhood = neighborhood,
                                            first = first,
                                            cdr3 = cdr3,
                                            progress = progress,
                                            verbose = verbose)
    
    # plot neighborhoods
    p2 <- list()
    p <- plotNeighborhoods(neighborhoods, 
                           threshold = threshold,
                           vline_size = vline_size,
                           binwidth=binwidth, center=center)
    p2[[length(p2)+1]] <- p
    
    # return all
    clonalAnalysisResult<-new("clonalAnalysisResult",
                              interVsIntera=interVsIntera,
                              threshold=threshold,
                              meanInter=meanInter,
                              sdInter=sdInter,
                              meanIntera=meanIntera,
                              sdIntera=sdIntera,
                              plotInterVsIntera=p1,
                              neighborhoods=neighborhoods,
                              plotNeighborhoods=p2)
    return(clonalAnalysisResult)        
}