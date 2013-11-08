setClass(
    "GeneOverlap", 
    representation(listA="character",
                   listB="character",
                   intersection="character",
                   union="character",
                   genome.size="numeric",
                   cont.tbl="matrix",
                   odds.ratio="numeric",
                   pval="numeric",
                   Jaccard="numeric",
                   is.tested="logical"),
    validity=function(object) {
        if(length(object@listA) > 0 && is.na(object@listA)) {
            stop("listA cannot be NA. Check your input.")
        }
        if(length(object@listB) > 0 && is.na(object@listB)) {
            stop("listB cannot be NA. Check your input.")
        }
        
        if(length(object@union) > object@genome.size) {
            stop("Union must NOT be larger than genome size.")
        }
    }
)

# Constructor
newGeneOverlap <- function(listA, listB, genome.size=NULL, 
                           spec=c('mm9.gene', 'hg19.gene', 'rn4.gene')) {
    listA <- unique(as.character(listA))
    listB <- unique(as.character(listB))
    listA <- listA[!is.na(listA)]
    listB <- listB[!is.na(listB)]
    
    # Setup genome size.
    if(is.null(genome.size)){
        spec <- match.arg(spec)
        genome.size <- switch(spec, 
                              mm9.gene=23000, 
                              hg19.gene=25000, 
                              rn4.gene=17000)
    }
    genome.size <- as.integer(genome.size)
    
    new("GeneOverlap", listA=listA, listB=listB, 
        intersection=intersect(listA, listB),
        union=union(listA, listB), 
        genome.size=genome.size, is.tested=F)
}

# Display methods.
setMethod(
    "show", "GeneOverlap",
    function(object) {
        cat("GeneOverlap object:\n")
        cat(sprintf("listA size=%d\n", length(object@listA)))
        cat(sprintf("listB size=%d\n", length(object@listB)))
        cat(sprintf("Intersection size=%d\n", length(object@intersection)))
        if(object@is.tested) {
            cat(sprintf("Overlapping p-value=%s\n", 
                        ifelse(object@pval < .01, 
                               format(object@pval, scientific=T, digits=2),
                               format(object@pval, digits=2)
                        )))
            cat(sprintf("Jaccard Index=%.1f\n", object@Jaccard))
        } else {
            cat("Overlap testing has not been performed yet.\n")
        }
    }
)

setMethod(
    "print", "GeneOverlap",
    function(x, ...) {
        cat("Detailed information about this GeneOverlap object:\n")
        cat(sprintf("listA size=%d, e.g. %s\n", 
                    length(x@listA), 
                    paste(head(x@listA), collapse=" ")))
        cat(sprintf("listB size=%d, e.g. %s\n", 
                    length(x@listB), 
                    paste(head(x@listB), collapse=" ")))
        cat(sprintf("Intersection size=%d, e.g. %s\n", 
                    length(x@intersection),
                    paste(head(x@intersection), collapse=" ")))
        cat(sprintf("Union size=%d, e.g. %s\n", 
                    length(x@union),
                    paste(head(x@union), collapse=" ")))
        cat(sprintf("Genome size=%d\n", x@genome.size))
        if(x@is.tested) {
            cat("# Contingency Table:\n")
            print(x@cont.tbl)
            cat(sprintf("Overlapping p-value=%s\n", 
                        ifelse(x@pval < .01, 
                               format(x@pval, scientific=T, digits=2),
                               format(x@pval, digits=2)
                        )))
            cat(sprintf("Odds ratio=%.1f\n", x@odds.ratio))
            cat("Overlap tested using Fisher's exact test (alternative=greater)\n")
            cat(sprintf("Jaccard Index=%.1f\n", x@Jaccard))
        } else {
            cat("Overlap has not been tested yet. Use testGeneOverlap method.\n")
        }
    }
)


# Accessors.
setGeneric("getListA", 
           function(object) { standardGeneric("getListA")}
           )
setMethod(
    "getListA", "GeneOverlap",
    function(object) {
        object@listA
    }
)

setGeneric("getListB", 
           function(object) { standardGeneric("getListB")}
)
setMethod(
    "getListB", "GeneOverlap",
    function(object) {
        object@listB
    }
)

setGeneric("getIntersection", 
           function(object) { standardGeneric("getIntersection")}
)
setMethod(
    "getIntersection", "GeneOverlap",
    function(object) {
        object@intersection
    }
)

setGeneric("getUnion", 
           function(object) { standardGeneric("getUnion")}
)
setMethod(
    "getUnion", "GeneOverlap",
    function(object) {
        object@union
    }
)

setGeneric("getGenomeSize", 
           function(object) { standardGeneric("getGenomeSize")}
)
setMethod(
    "getGenomeSize", "GeneOverlap",
    function(object) {
        object@genome.size
    }
)

setGeneric("getTested", 
           function(object) { standardGeneric("getTested")}
)
setMethod(
    "getTested", "GeneOverlap",
    function(object) {
        object@is.tested
    }
)

setGeneric("getContbl", 
           function(object) { standardGeneric("getContbl")}
)
setMethod(
    "getContbl", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@cont.tbl
        } else {
            warning("Test has not been performed yet.\n")
            matrix(nrow=0, ncol=0)
        }
    }
)

setGeneric("getPval", 
           function(object) { standardGeneric("getPval")}
)
setMethod(
    "getPval", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@pval
        } else {
            warning("Test has not been performed yet.\n")
            NA
        }
    }
)

setGeneric("getOddsRatio", 
           function(object) { standardGeneric("getOddsRatio")}
)
setMethod(
    "getOddsRatio", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@odds.ratio
        } else {
            warning("Test has not been performed yet.\n")
            NA
        }
    }
)

setGeneric("getJaccard", 
           function(object) { standardGeneric("getJaccard")}
)
setMethod(
    "getJaccard", "GeneOverlap",
    function(object) {
        if(object@is.tested) {
            object@Jaccard
        } else {
            warning("Test has not been performed yet.\n")
            NA
        }
    }
)

setGeneric("setListA<-", 
           function(object, value) { standardGeneric("setListA<-") }
)
setReplaceMethod(
    "setListA", "GeneOverlap",
    function(object, value) {
        object@listA <- as.character(value)
        object@intersection <- intersect(object@listA, object@listB)
        object@union <- union(object@listA, object@listB)
        object@is.tested <- F
        validObject(object)
        
        object
    }
)

setGeneric("setListB<-", 
           function(object, value) { standardGeneric("setListB<-") }
)
setReplaceMethod(
    "setListB", "GeneOverlap",
    function(object, value) {
        object@listB <- as.character(value)
        object@intersection <- intersect(object@listA, object@listB)
        object@union <- union(object@listA, object@listB)
        object@is.tested <- F
        validObject(object)
        
        object
    }
)

setGeneric("setGenomeSize<-", 
           function(object, value) { standardGeneric("setGenomeSize<-") }
)
setReplaceMethod(
    "setGenomeSize", "GeneOverlap",
    function(object, value) {
        object@genome.size <- value
        object@is.tested <- F
        validObject(object)
        
        object
    }
)


# Test function.
setGeneric("testGeneOverlap", 
           function(object) { standardGeneric("testGeneOverlap") }
           )
setMethod(
    "testGeneOverlap", "GeneOverlap",
    function(object) {
        # Configure contingency table.
        sizeA <- length(object@listA)
        sizeB <- length(object@listB)
        object@cont.tbl <- matrix(c(object@genome.size - length(object@union), 
                                    sizeB - length(object@intersection), 
                                    sizeA - length(object@intersection), 
                                    length(object@intersection)), 
                                  ncol=2)
        rownames(object@cont.tbl) <- c('notB', 'inB')
        colnames(object@cont.tbl) <- c('notA', 'inA')
        
        # Perform Fisher's exact test.
        res.fisher <- try(fisher.test(object@cont.tbl, alternative='greater'), 
                          silent=TRUE)
        if(is.list(res.fisher)) {
            object@odds.ratio <- setNames(res.fisher$estimate, NULL)
            object@pval <- res.fisher$p.value
        } else {
            object@odds.ratio <- .0
            object@pval <- 1.
        }
        
        # Calculate Jaccard index.
        object@Jaccard <- ifelse(length(object@union) == 0, 0, 
                                 length(object@intersection) / 
                                     length(object@union)
                                 )
        
        object@is.tested <- T
        object
    }
)

