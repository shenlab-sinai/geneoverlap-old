\name{testGeneOverlap}
\alias{testGeneOverlap-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test function for the GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
Perform Fisher's exact test based on the information supplied in the GeneOverlap object, such as gene list A, B and genome size. Set the tested Boolean label.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
testGeneOverlap(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A GeneOverlap object.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A GeneOverlap object with valid p-value, odds ratio and contingency table. The tested Boolean label is set to true. Use show or print to display a summary of the object. Use accessors to get information of each slot.
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
getOddsRatio(go.obj)
getContbl(go.obj)
print(go.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
