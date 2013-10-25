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



ExtractAsMatrix <- function(m, name) {
# Extract information from a GeneOverlapMatrix object by designated name as 
# a matrix.
    
    sapply(m, function(ci) {
        sapply(ci, function(ri) {
            ri[[name]]
        })
    })
}

print.GeneOverlapMatrix <- function(x, ...) {
    or.mat <- ExtractAsMatrix(x$overlap.matrix, "odds.ratio")
    pv.mat <- ExtractAsMatrix(x$overlap.matrix, "pval")
    it.mat <- ExtractAsMatrix(x$overlap.matrix, "intersection")
 
    print(list(odds.ratio=or.mat, pval=pv.mat, intersection=it.mat))
}

summary.GeneOverlapMatrix <- function(obj, ...) {
    summ.res <- list()
    summ.res$desc <- paste("Overlap matrix based on: ", 
                           ifelse(obj$self.compare, 
                                  "self-comparison on one gene set.",
                                  "cross-comparison on two gene sets."),
                           sep="")
    summ.res$genome.size <- obj$overlap.matrix[[1]][[1]]$genome.size
    
    or.v <- c(ExtractAsMatrix(obj$overlap.matrix, "odds.ratio"))
    summ.res$or.v <- or.v[or.v > 1]
    
    pv.v <- c(ExtractAsMatrix(obj$overlap.matrix, "pval"))
    summ.res$pv.v <- pv.v[pv.v < 1]
    
    it.v <- c(ExtractAsMatrix(obj$overlap.matrix, "intersection"))
    summ.res$it.v <- it.v[it.v > 0]
    
    structure(summ.res, class="summary.GeneOverlapMatrix")
}

print.summary.GeneOverlapMatrix <- function(x, ...) {
    cat(sprintf("Description -\n%s\n\n", x$desc))
    cat(sprintf("Genome size: %d\n\n", x$genome.size))
    cat(sprintf("%d odds ratios are larger than 1:\n", length(x$or.v)))
    print(summary(x$or.v))
    cat(sprintf("\n%d p-values are less than 1:\n", length(x$pv.v)))
    print(summary(x$pv.v))
    cat(sprintf("\n%d intersections are larger than 0:\n", length(x$it.v)))
    print(summary(x$it.v))
}

plot.GeneOverlapMatrix <- function(x, adj.p=F, cutoff=.05, ncolused=9, 
                                   grid.col=c("Greens", "Blues", "Greys", 
                                              "Oranges", "Purples", "Reds"),
                                   note.col="red", ...) {
    
    or.mat <- ExtractAsMatrix(x$overlap.matrix, "odds.ratio")
    pv.mat <- ExtractAsMatrix(x$overlap.matrix, "pval")

    # Adjust p-values if needed.
    if(adj.p) {
        n <- ifelse(x$self.compare, 
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
    grid.col <- match.arg(grid.col)
    
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














