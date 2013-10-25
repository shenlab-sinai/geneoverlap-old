\name{GeneOverlapMatrix}
\alias{GeneOverlapMatrix-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Represent gene list overlaps as a GeneOverlap matrix
}
\description{
Given one or two gene sets of gene lists, create a GeneOverlap matrix to represent all pairwise comparisons. This S3 class also provides functions to print, summary and plot the matrix.
}
\usage{
gom <- GeneOverlapMatrix(gsetA, gsetB = NULL, ...)
print(gom)
summary(gom)
plot(gom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsetA}{
  Gene set A as a named list with each element being a vector of gene names.
}
  \item{gsetB}{
  Gene set B as a named list with each element being a vector of gene names.
}
  \item{\dots}{
  These arguments will be passed on to the GeneOverlap constructor.
}
  \item{gom}{
  A GeneOverlapMatrix object.
}
  \item{adj.p}{
  Boolean tag for whether p-values should be ajusted.
}
  \item{cutoff}{
  P-value cutoff to mask the insignificant comparisons.
}
  \item{ncolused}{
  Number of colors used to represent the scale of odds ratios.
}
  \item{grid.col}{
  Color for odds ratios.
}
  \item{note.col}{
  Color for p-value text.
}
}
\details{
This will create a matrix so that each grid represents the pairwise comparison between two gene lists. Given two sets A and B, the matrix will represent all comparisons bewteen gene lists in A vs. gene lists in B. The set A will be shown as rows and the set B will be shown as columns. If given only one set, it will represent all comparisons between gene lists within the set and only the upper triangular of the matrix will be used.

The grids that are below the p-value cutoff or have odds ratio no larger than 1 will be shown the lightest color.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of lists of GeneOverlap objects. The outer list represents the columns of the matrix and the inner list represents the rows of the matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
\url{http://en.wikipedia.org/wiki/Fisher's_exact_test}

%% ~put references to the literature/web site here ~
}
\author{
Li Shen <\email{shenli.sam@gmail.com}>
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GeneOverlap-class}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Use random numbers to create some fake gene lists.
gsetA <- setNames(lapply(1:2, function(x) 
        as.character(sample(1:100, size=20))), c("A1", "A2"))
gsetB <- setNames(lapply(1:3, function(x) 
        as.character(sample(1:100, size=20))), c("B1", "B2", "B3"))
gom <- GeneOverlapMatrix(gsetA, gsetB, genome.size=100)
print(gom)  # odds ratios, p-values and intersections.
summary(gom)
## Make sure the plotting area is large enough to run the following code!
# plot(gom, cutoff=.1)  # represent the overlaps as a heatmap.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ graphs }% __ONLY ONE__ keyword per line
