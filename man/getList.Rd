\name{getList}
\alias{getListA}
\alias{getListA-method}
\alias{setListA<-}
\alias{setListA<--method}
\alias{getListB}
\alias{getListB-method}
\alias{setListB<-}
\alias{setListB<--method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessors for the "listA" and "listB" slots of GeneOverlap class
%%  ~~function to do ... ~~
}
\description{
The listA and listB slots hold the gene lists A and B as character vectors.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getListA(object)
setListA(object) <- gene.list.A
getListB(object)
setListB(object) <- gene.list.B
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A GeneOverlap object.
%%     ~~Describe \code{x} here~~
}
}
\details{
After setListX function is called, the tested Boolean label will be reset to false.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{listA }{A character vector represents gene list A.}
\item{listB }{A character vector represents gene list B.}
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
\code{\link{GeneOverlap}}, \code{\link{newGeneOverlap}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
listA <- c("A", "B", "C")
listB <- c("B", "C", "D")
go.obj <- newGeneOverlap(listA, listB, genome.size=10)
go.obj <- testGeneOverlap(go.obj)
getListB(go.obj)
getTested(go.obj)  # true.
setListB(go.obj) <- c("C", "D", "E")
getTested(go.obj) # false.
go.obj <- testGeneOverlap(go.obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
