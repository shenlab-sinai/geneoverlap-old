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
                    if(ri == ci) {
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
           function(object, name) { standardGeneric("getMatrix")}
)
setMethod(
    "getMatrix", "GeneOverlapMatrix",
    function(object, name=c("pval", "odds.ratio", "intersection", "union")) {
        name <- match.arg(name)
        sapply(object@go.nested.list, function(ci) {
            sapply(ci, function(ri) {
                switch(name, 
                       pval=getPval(ri), 
                       odds.ratio=getOddsRatio(ri),
                       intersection=length(getIntersection(ri)),
                       union=length(getUnion(ri)))
            })
        })
    }
)

setGeneric("getNestedList", 
           function(object, name) { standardGeneric("getNestedList")}
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
            stopifnot(j >= 1 && j <= length(x@go.nested.list))
        } else {
            stopifnot(j %in% names(x@go.nested.list))
        }
        gom.col <- x@go.nested.list[[j]]
        if(is.numeric(i)) {
            i <- as.integer(i)
            stopifnot(i >= 1 && i <= length(gom.col))
        } else {
            stopifnot(i %in% names(gom.col))
        }
        gom.col[[i]]
    }
)

# Visualization function.
setGeneric("drawHeatmap", 
           function(object, adj.p=F, cutoff=.05, ncolused=9, 
                    grid.col=c("Greens", "Blues", "Greys", 
                               "Oranges", "Purples", "Reds"),
                    note.col="red") { 
               standardGeneric("drawHeatmap") 
           }
)
setMethod(
    "drawHeatmap", "GeneOverlapMatrix",
    function(object, adj.p=F, cutoff=.05, ncolused=9, 
             grid.col=c("Greens", "Blues", "Greys", 
                        "Oranges", "Purples", "Reds"),
             note.col="red") {
        
        grid.col <- match.arg(grid.col)
        or.mat <- getMatrix(object, "odds.ratio")
        pv.mat <- getMatrix(object, "pval")
        
        # Adjust p-values if needed.
        if(adj.p) {
            n <- ifelse(object@self.compare, 
                        (nrow(pv.mat) + 1) * nrow(pv.mat) / 2,
                        nrow(pv.mat) * ncol(pv.mat))
            pv.mat <- matrix(p.adjust(pv.mat, method='BH', n=n), 
                             nrow=nrow(pv.mat))
        }
        
        # Use odds ratio and p-value cutoff to mask insignificant cells.
        plot.mat <- or.mat
        plot.mat[ pv.mat >= cutoff | or.mat < 1.0 ] <- 1.0
        
        # Cell notes.
        note.mat <- format(pv.mat, digits=1)
        note.mat[pv.mat < .01] <- format(pv.mat, digits=1, 
                                         scientific=T)[pv.mat < .01]
        note.mat[plot.mat == 1] <- ''
        
        # Configure heatmap graphic properties.
        row_sep <- 1:(nrow(plot.mat) - 1)
        col_sep <- 1:(ncol(plot.mat) - 1)
        longedge <- max(nrow(plot.mat), ncol(plot.mat))
        row_cexrc <- 0.4 + 1/log10(longedge + 2)
        col_cexrc <- row_cexrc
        key_size <- 0.2 + 1 / log10(longedge + 4)
        margins_use <- c(max(nchar(colnames(plot.mat))) * 0.8 + 5, 
                         max(nchar(rownames(plot.mat))) * 0.8 + 5)
        
        # Draw the heatmap!
        heatmap.2(plot.mat, cellnote=note.mat, 
                  col=brewer.pal(ncolused, grid.col), notecol=note.col, 
                  margins=margins_use, colsep=col_sep, rowsep=row_sep, 
                  key=T, keysize=key_size,
                  cexRow=row_cexrc, cexCol=col_cexrc, 
                  scale='none', Colv=NA, Rowv=NA, trace='none', 
                  dendrogram='none', density.info='none', 
                  sepcolor='white', sepwidth=c(0.002,0.002),
                  notecex=1.8)
    }
)

















