\name{getGenomeSize}
\alias{getGenomeSize-method}
\alias{setGenomeSize<-}
\alias{setGenomeSize<--method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessors for the "genome.size" slot of the GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
The genome.size slot contains the number of genes in the genome as an integer.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getGenomeSize(object)
setGenomeSize(object) <- new.genome.size
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A GeneOverlap object.
%%     ~~Describe \code{x} here~~
}
}
\details{
After setGenomeSize function is called, the tested Boolean label will be reset to false.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{genome.size }{An integer represents the genome size.}
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
getPval(go.obj)
setGenomeSize(go.obj) <- 20
getPval(go.obj)  # NA.
# See how changing the genome size will affect the p-value.
v.gs <- c(10, 100, 1000, 10e3, 100e3)
v.gs <- setNames(v.gs, as.character(v.gs))
sapply(v.gs, function(x) {
    setGenomeSize(go.obj) <- x
    getPval(testGeneOverlap(go.obj))
})
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
