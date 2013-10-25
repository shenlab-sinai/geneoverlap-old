\name{getReadonly}
\alias{getIntersection}
\alias{getIntersection-method}
\alias{getUnion}
\alias{getUnion-method}
\alias{getTested}
\alias{getTested-method}
\alias{getContbl}
\alias{getContbl-method}
\alias{getPval}
\alias{getPval-method}
\alias{getOddsRatio}
\alias{getOddsRatio-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read-only accessors for the "intersection", "union", "is.tested", "cont.tbl", "pval", "odds.ratio" slots of the GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
The intersection and union slots contain the gene names as character vectors. The is.tested slot contains a Boolean label indicating whether the object has been tested or not. The cont.tbl slot contains the contingency table as a matrix. The pval and odds.ratio slots contain the p-value and odds ratio as numerics, respectively. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getIntersection(object)
getUnion(object)
getTested(object)
getContbl(object)
getPval(object)
getOddsRatio(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A GeneOverlap object.
%%     ~~Describe \code{x} here~~
}
}
\details{
If the GeneOverlap object has not been tested yet, the returned p-value and odds ratio will be NA, the contingency table will be an empty matrix.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{intersection }{A character vector represents the overlapped genes.}
\item{union }{A character vector represents the genes in the union of A and B.}
\item{is.tested }{A Boolean represents whether the overlapping test has been performed or not.}
\item{cont.tbl }{A matrix represents the contingency table.}
\item{pval }{A numeric represents the significance of overlap.}
\item{odds.ratio }{A numeric represents the odds ratio in comparison to the genomic background.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
listA <- c("A", "B", "C")
listB <- c("B", "C", "D")
go.obj <- newGeneOverlap(listA, listB, genome.size=10)
go.obj <- testGeneOverlap(go.obj)
getIntersection(go.obj)
getUnion(go.obj)
getTested(go.obj)
getContbl(go.obj)
getPval(go.obj)
getOddsRatio(go.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
