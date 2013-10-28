\name{newGOM}
\alias{newGOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructor for the GeneOverlapMatrix class
%%  ~~function to do ... ~~
}
\description{
Use this function to create objects of the GeneOverlapMatrix class.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
newGOM(gsetA, gsetB=list(), genome.size=NULL, 
    spec=c('mm9.gene', 'hg19.gene', 'rn4.gene'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsetA}{
  Gene set A as a named list with each element being a vector of gene names.
}
  \item{gsetB}{
  Gene set B as a named list with each element being a vector of gene names.
}
  \item{genome.size}{
  The number of genes in the genome as an integer.
}
  \item{spec}{
  A string description of the genome with a few presetted genome sizes.
}
%%     ~~Describe \code{x} here~~
}
\details{
This will create a matrix so that each grid represents the pairwise comparison between two gene lists. Given two sets A and B, the matrix will represent all comparisons bewteen gene lists in A vs. gene lists in B. The set A will be shown as rows and the set B will be shown as columns. If given only one set, it will represent all comparisons between gene lists within the set and only the upper triangular of the matrix will be used.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
