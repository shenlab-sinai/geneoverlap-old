\name{GeneOverlapMatrix}
\alias{GeneOverlapMatrix-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix representation of the pairwise overlaps between two gene sets
}
\description{
Given one or two gene sets each contains one or more gene lists, create a matrix to represent all pairwise comparisons between the gene lists. This class also provides functions to visualize the matrix.
}
\details{
The problem is stated as the representation of all pairwise comparisons between two gene sets each contains one or more gene lists. This is represented as a matrix where the rows correspond to one gene set and the columns correspond to the other. Each grid represents the overlap between two gene lists of the corresponding row and column. This class calls the GeneOverlap constructor to create objects that represent the overlapping information. When there is only one gene set, the matrix represents the self-comparison within the gene set and only the upper triangle is used.

The significance of gene overlap is characterized by two pieces of information: odds ratio and p-value. This class provides functions to visualize these information as a heatmap. The darkness of the color of each grid represents the odds ratio while the texts superimposed on the grids state the p-values.
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
Li Shen <\email{shenli.sam@gmail.com}>

\url{http://www.linkedin.com/in/lshen/}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ hplot }
\keyword{ graphs }% __ONLY ONE__ keyword per line
