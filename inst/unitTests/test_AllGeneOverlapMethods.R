test_GeneOverlap <- function() {
    # Vanilla test.
    listA <- c("A", "B", "A")
    listB <- c("B", "C", "B")
    manual.contbl <- matrix(c(7, 1, 1, 1), nrow=2)
    fish.res <- fisher.test(manual.contbl, alternative="greater")
    go.obj <- newGeneOverlap(listA, listB, genome.size=10)
    checkEquals(getListA(go.obj), c("A", "B"))
    checkEquals(getListB(go.obj), c("B", "C"))
    checkEqualsNumeric(getGenomeSize(go.obj), 10)
    checkEqualsNumeric(length(getIntersection(go.obj)), 1)
    checkEqualsNumeric(length(getUnion(go.obj)), 3)
    checkEqualsNumeric(getPval(go.obj), NA)
    checkEqualsNumeric(getOddsRatio(go.obj), NA)
    go.obj <- testGeneOverlap(go.obj)
    checkEqualsNumeric(getPval(go.obj), fish.res$p.value)
    checkEqualsNumeric(getOddsRatio(go.obj), fish.res$estimate)
    
    # Test NO overlap.
    listA <- c("A", "B")
    listB <- c("C", "D")
    go.obj <- testGeneOverlap(newGeneOverlap(listA, listB, genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    
    # Test absolute overlap (p-value=0).
    listA <- LETTERS
    listB <- LETTERS
    go.obj <- testGeneOverlap(newGeneOverlap(listA, listB, genome.size=100))
    checkEqualsNumeric(getPval(go.obj), 0)
    checkEqualsNumeric(getOddsRatio(go.obj), Inf)
    
    # Test empty gene list.
    go.obj <- testGeneOverlap(newGeneOverlap("A", NULL, genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    go.obj <- testGeneOverlap(newGeneOverlap(NULL, "B", genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    go.obj <- testGeneOverlap(newGeneOverlap(NULL, NULL, genome.size=10))
    checkEqualsNumeric(getPval(go.obj), 1)
    checkEqualsNumeric(getOddsRatio(go.obj), 0)
    
    # Unknown species.
    checkException(newGeneOverlap("A", "B", spec="speciesdoesnotexist"))

    # Genome smaller than gene lists combined.
    checkException(newGeneOverlap("A", "B", genome.size=1))

    # Bad gene lists.
    checkException(newGeneOverlap("A", NA))
    checkException(newGeneOverlap(NA, "B"))
}

test_GeneOverlapMatrix <- function() {
    # Vanilla test.
    gv1 <- c("A", "B")
    gv2 <- c("B", "C")
    manual.contbl <- matrix(c(7, 1, 1, 1), nrow=2)
    fish.res <- fisher.test(manual.contbl, alternative="greater")
    gsetA <- list(A=gv1)
    gsetB <- list(B=gv2)
    gom.obj <- newGOM(gsetA, gsetB, genome.size=10)
    checkEquals(getGsetA(gom.obj), gsetA)
    checkEquals(getGsetB(gom.obj), gsetB)
    checkEquals(getSelfCompare(gom.obj), F)
    checkEqualsNumeric(getMatrix(gom.obj, "pval"), fish.res$p.value)
    checkEqualsNumeric(getMatrix(gom.obj, "odds.ratio"), fish.res$estimate)
    checkEqualsNumeric(getMatrix(gom.obj, "intersection"), 1)
    checkEqualsNumeric(getMatrix(gom.obj, "union"), 3)
    checkEquals(getNestedList(gom.obj, "intersection")[[1]][[1]], "B")
    checkEquals(getNestedList(gom.obj, "union")[[1]][[1]], c("A", "B", "C"))
    go.obj <- testGeneOverlap(newGeneOverlap(gv1, gv2, genome.size=10))
    checkEquals(gom.obj[1, 1], go.obj)
    checkEquals(gom.obj["A", "B"], go.obj)
    
}

#     # gsetA not long enough.
#     gs1 <- list(A=c("A"))
#     checkException(GeneOverlapMatrix(gs1))
#     # gsetB is empty.
#     gs2 <- list()
#     checkException(GeneOverlapMatrix(gs1, gs2))
#     
#     # A simple case of two sets of one gene list each with no overlap.
#     gs3 <- list(B=c("B"))
#     obj.no.overlap <- structure(list(pval=1.0, odds.ratio=.0, intersection=0, 
#                                      sizeA=1, sizeB=1, union=2, genome.size=20e3),
#                                 class="GeneOverlap")
#     mat.no.overlap <- structure(list(self.compare=F, 
#                                      overlap.matrix=list(
#                                          B=list(A=obj.no.overlap))
#                                      ),
#                                 class="GeneOverlapMatrix")
#     checkEquals(GeneOverlapMatrix(gs1, gs3, genome.size=20e3), mat.no.overlap)
    










