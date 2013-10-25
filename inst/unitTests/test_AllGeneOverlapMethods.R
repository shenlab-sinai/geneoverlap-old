test_GeneOverlap <- function() {
    # Vanilla test.
    listA <- c("A", "B", "A")
    listB <- c("B", "C", "B")
    go.obj <- newGeneOverlap(listA, listB, genome.size=10)
    checkEquals(getListA(go.obj), c("A", "B"))
    checkEquals(getListB(go.obj), c("B", "C"))
    checkEquals(getGenomeSize(go.obj), 10)
    checkEquals(length(getIntersection(go.obj)), 1)
    checkEquals(length(getUnion(go.obj)), 3)
    checkEquals(getPval(go.obj), NA)
    checkEquals(getOddsRatio(go.obj), NA)
    go.obj <- testGeneOverlap(go.obj)
    checkEquals(getPval(go.obj), 0.3777778, tolerance=1e-5)
    checkEquals(getOddsRatio(go.obj), 5.291552, tolerance=1e-5)
    
    # Test NO overlap.
    listA <- c("A", "B")
    listB <- c("C", "D")
    go.obj <- testGeneOverlap(newGeneOverlap(listA, listB, genome.size=10))
    checkEquals(getPval(go.obj), 1)
    checkEquals(getOddsRatio(go.obj), 0)
    
    # Test absolute overlap (p-value=0).
    listA <- LETTERS
    listB <- LETTERS
    go.obj <- testGeneOverlap(newGeneOverlap(listA, listB, genome.size=100))
    checkEquals(getPval(go.obj), 0)
    checkEquals(getOddsRatio(go.obj), Inf)
    
    # Test empty gene list.
    go.obj <- testGeneOverlap(newGeneOverlap("A", NULL, genome.size=10))
    checkEquals(getPval(go.obj), 1)
    checkEquals(getOddsRatio(go.obj), 0)
    go.obj <- testGeneOverlap(newGeneOverlap(NULL, "B", genome.size=10))
    checkEquals(getPval(go.obj), 1)
    checkEquals(getOddsRatio(go.obj), 0)
    go.obj <- testGeneOverlap(newGeneOverlap(NULL, NULL, genome.size=10))
    checkEquals(getPval(go.obj), 1)
    checkEquals(getOddsRatio(go.obj), 0)
    
    # Unknown species.
    checkException(newGeneOverlap("A", "B", spec="speciesdoesnotexist"))

    # Genome smaller than gene lists combined.
    checkException(newGeneOverlap("A", "B", genome.size=1))

    # Bad gene lists.
    checkException(newGeneOverlap("A", NA))
    checkException(newGeneOverlap(NA, "B"))
}

test_GeneOverlapMatrix <- function() {
    # gsetA not long enough.
    gs1 <- list(A=c("A"))
    checkException(GeneOverlapMatrix(gs1))
    # gsetB is empty.
    gs2 <- list()
    checkException(GeneOverlapMatrix(gs1, gs2))
    
    # A simple case of two sets of one gene list each with no overlap.
    gs3 <- list(B=c("B"))
    obj.no.overlap <- structure(list(pval=1.0, odds.ratio=.0, intersection=0, 
                                     sizeA=1, sizeB=1, union=2, genome.size=20e3),
                                class="GeneOverlap")
    mat.no.overlap <- structure(list(self.compare=F, 
                                     overlap.matrix=list(
                                         B=list(A=obj.no.overlap))
                                     ),
                                class="GeneOverlapMatrix")
    checkEquals(GeneOverlapMatrix(gs1, gs3, genome.size=20e3), mat.no.overlap)
    
}









