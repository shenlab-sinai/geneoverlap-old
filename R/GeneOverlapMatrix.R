setClass(
    "GeneOverlapMatrix", 
    representation(gsetA="list",
                   gsetB="list",
                   self.compare="logical",
                   go.nested.list="list"),
    validity=function(object) {
        if(length(object@gsetB) == 0) {
            stopifnot(length(object@gsetA) > 1 && object@self.compare)
        } else {
            stopifnot(length(object@gsetA) > 0 && !object@self.compare)
        }

    }
)

# Constructor.
newGOM <- function(gsetA, gsetB=list(), genome.size=NULL, 
                   spec=c('mm9.gene', 'hg19.gene', 'rn4.gene')) {
    stopifnot(is.list(gsetA) && is.list(gsetB))
    # Construct GeneOverlap objects for all pairwise comparisons.
    if(length(gsetB) == 0) {
        stopifnot(length(gsetA) >= 2)
        self.compare <- T
        row.iter <- 1:(length(gsetA) - 1)
        col.iter <- 2:length(gsetA)
        go.nested.list <- 
            lapply(col.iter, function(ci) {
                this.col <- lapply(row.iter, function(ri) {
                    if(ri >= ci) {
                        go.obj <- newGeneOverlap(NULL, NULL)  # same list.
                        testGeneOverlap(go.obj)
                    } else {
                        go.obj <- newGeneOverlap(gsetA[[ri]], gsetA[[ci]], 
                                                 genome.size, spec)
                        testGeneOverlap(go.obj)
                    }
                })
                names(this.col) <- names(gsetA)[row.iter]
                this.col
        })
        names(go.nested.list) <- names(gsetA)[col.iter]
    } else {
        stopifnot(length(gsetA) >= 1 && length(gsetB) >= 1)
        self.compare <- F
        go.nested.list <- 
            lapply(gsetB, function(b) {
                this.col <- lapply(gsetA, function(a) {
                    go.obj <- newGeneOverlap(a, b, genome.size, spec)
                    testGeneOverlap(go.obj)
                })
                names(this.col) <- names(gsetA)
                this.col
            })
        names(go.nested.list) <- names(gsetB)
    }
    
    new("GeneOverlapMatrix", gsetA=gsetA, gsetB=gsetB, 
        self.compare=self.compare, go.nested.list=go.nested.list)
}

# show and print functions.
setMethod("show", "GeneOverlapMatrix",
          function(object) {
              gom.dim <- dim(getMatrix(object, "pval"))
              cat(sprintf("A <%d x %d> GeneOverlapMatrix object\n", 
                          gom.dim[1], gom.dim[2]))
              gsetA <- getGsetA(object)
              gsetB <- getGsetB(object)
              cat("Geneset A sizes:\n")
              print(sapply(gsetA, length))
              if(getSelfCompare(object)) {
                  cat("Matrix is based on self-comparison of geneset A.\n")
              } else {
                  cat("Geneset B sizes:\n")
                  print(sapply(gsetB, length))
              }
          }
)

setMethod("print", "GeneOverlapMatrix",
          function(x, ...) {
              cat("A GeneOverlapMatrix object:\n")
              
              int.mat <- getMatrix(x, "intersection")
              cat("###### Intersection ######\n")
              print(int.mat)
              
              pval.mat <- getMatrix(x, "pval")
              cat("###### P-value ######\n")
              print(pval.mat)
              
              or.mat <- getMatrix(x, "odds.ratio")
              cat("###### Odds Ratio ######\n")
              print(or.mat)
              
              ja.mat <- getMatrix(x, "Jaccard")
              cat("###### Jaccard Index ######\n")
              print(ja.mat)
          }
)


# Accessors.
setGeneric("getGsetA", 
           function(object) { standardGeneric("getGsetA")}
)
setMethod(
    "getGsetA", "GeneOverlapMatrix",
    function(object) {
        object@gsetA
    }
)

setGeneric("getGsetB", 
           function(object) { standardGeneric("getGsetB")}
)
setMethod(
    "getGsetB", "GeneOverlapMatrix",
    function(object) {
        object@gsetB
    }
)

setGeneric("getSelfCompare", 
           function(object) { standardGeneric("getSelfCompare")}
)
setMethod(
    "getSelfCompare", "GeneOverlapMatrix",
    function(object) {
        object@self.compare
    }
)

setGeneric("getMatrix", 
           function(object, name=c("pval", "odds.ratio", "intersection", 
                                   "union", "Jaccard")) { 
               standardGeneric("getMatrix")
           }
)
setMethod(
    "getMatrix", "GeneOverlapMatrix",
    function(object, name=c("pval", "odds.ratio", "intersection", "union", 
                            "Jaccard")) {
        name <- match.arg(name)
        sapply(object@go.nested.list, function(ci) {
            sapply(ci, function(ri) {
                switch(name, 
                       pval=getPval(ri), 
                       odds.ratio=getOddsRatio(ri),
                       intersection=length(getIntersection(ri)),
                       union=length(getUnion(ri)),
                       Jaccard=getJaccard(ri)
                       )
            })
        })
    }
)

setGeneric("getNestedList", 
           function(object, name=c("intersection", "union", 
                                   "cont.tbl")) { 
               standardGeneric("getNestedList")
           }
)
setMethod(
    "getNestedList", "GeneOverlapMatrix",
    function(object, name=c("intersection", "union", "cont.tbl")) {
        name <- match.arg(name)
        lapply(object@go.nested.list, function(ci) {
            lapply(ci, function(ri) {
                switch(name, 
                       intersection=getIntersection(ri),
                       union=getUnion(ri),
                       cont.tbl=getContbl(ri))
            })
        })
    }
)

# More complex accessors.
setMethod(
    "[", "GeneOverlapMatrix",
    function(x, i, j) {
        stopifnot(is.numeric(i) || is.character(i))
        stopifnot(is.numeric(j) || is.character(j))
        if(is.numeric(j)) {
            j <- as.integer(j)
            stopifnot(abs(j) <= length(x@go.nested.list))
        } else {
            stopifnot(j %in% names(x@go.nested.list))
        }
        gom.col <- x@go.nested.list[[j]]
        if(is.numeric(i)) {
            i <- as.integer(i)
            stopifnot(abs(i) <= length(gom.col))
        } else {
            stopifnot(i %in% names(gom.col))
        }
        gom.col[[i]]
    }
)

# Visualization function.
setGeneric("drawHeatmap", 
           function(object, what=c("odds.ratio", "Jaccard"), log.scale=F, 
                    adj.p=F, cutoff=.05, ncolused=9, 
                    grid.col=c("Greens", "Blues", "Greys", 
                               "Oranges", "Purples", "Reds"),
                    note.col="red") { 
               standardGeneric("drawHeatmap") 
           }
)
setMethod(
    "drawHeatmap", "GeneOverlapMatrix",
    function(object, what=c("odds.ratio", "Jaccard"), log.scale=F, adj.p=F, 
             cutoff=.05, ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
                                                "Oranges", "Purples", "Reds"),
             note.col="red") {

        # Arguments setup.
        stopifnot(cutoff > 0 && cutoff <= 1)
        what <- match.arg(what)
        grid.col <- match.arg(grid.col)

        # Matrix values.
        pv.mat <- getMatrix(object, "pval")
        plot.mat <- switch(what, 
                           odds.ratio=getMatrix(object, "odds.ratio"),
                           Jaccard=getMatrix(object, "Jaccard")
                           )
        if(what == "odds.ratio" && log.scale) {
            plot.mat <- log2(plot.mat)
        }
        
        # Adjust p-values if needed.
        pv.mask <- NULL
        if(object@self.compare) {
            pv.mask <- sapply(1:ncol(pv.mat), function(j) {
                c(rep(T, j), rep(F, nrow(pv.mat) - j))
            })
        }
        if(adj.p) {
            if(object@self.compare) {
                pv.mat[pv.mask] <- p.adjust(pv.mat[pv.mask], method='BH')
            } else {
                pv.mat <- matrix(p.adjust(pv.mat, method='BH'), 
                                 nrow=nrow(pv.mat))
            }
        }
        
        # Marker value of insignificant events.
        insig.val <- 1
        if(what == "odds.ratio" && log.scale || what == "Jaccard") {
            insig.val <- 0
        }
        
        # Use p-value cutoff to mask insignificant cells.
        plot.mat[ pv.mat >= cutoff ] <- insig.val
        
        # Cell notes.
        note.mat <- format(pv.mat, digits=1)
        note.mat[pv.mat < .01] <- format(pv.mat, digits=1, 
                                         scientific=T)[pv.mat < .01]
        note.mat[plot.mat == insig.val] <- "N.S."
        if(object@self.compare) { note.mat[ !pv.mask ] <- "--" }
        
        # Configure heatmap graphic properties.
        row_sep <- 1:(nrow(plot.mat) - 1)
        col_sep <- 1:(ncol(plot.mat) - 1)
        longedge <- max(nrow(plot.mat), ncol(plot.mat))
        row_cexrc <- 0.4 + 1/log10(longedge + 2)
        col_cexrc <- row_cexrc
        key_size <- 0.2 + 1 / log10(longedge + 4)
        margins_use <- c(max(nchar(colnames(plot.mat))) * 0.8 + 5, 
                         max(nchar(rownames(plot.mat))) * 0.8 + 5)
        main.txt <- switch(what, 
                           odds.ratio=ifelse(log.scale, "log2(Odds Ratio)", 
                                             "Odds Ratio"), 
                           Jaccard="Jaccard Index")
        footnote <- "N.S.: Not Significant; --: Ignored"
        # sidenote <- sprintf("Log Scale=%s", log.scale)
        
        # Draw the heatmap!
        heatmap.2(plot.mat, cellnote=note.mat, 
                  main=main.txt, xlab=footnote, # ylab=sidenote,
                  col=brewer.pal(ncolused, grid.col), notecol=note.col, 
                  margins=margins_use, colsep=col_sep, rowsep=row_sep, 
                  key=T, keysize=key_size,
                  cexRow=row_cexrc, cexCol=col_cexrc, 
                  scale='none', Colv=NA, Rowv=NA, trace='none', 
                  dendrogram='none', density.info='none', 
                  sepcolor='white', sepwidth=c(0.002,0.002),
                  notecex=1.6)
    }
)

















