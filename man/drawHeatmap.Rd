\name{drawHeatmap}
\alias{drawHeatmap}
\alias{drawHeatmap,GeneOverlapMatrix-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Visualize GeneOverlapMatrix objects as heatmaps
%%  ~~function to do ... ~~
}
\description{
    Visualization function for GeneOverlapMatrix objects. Use the darkness of 
    the colors to represent the odds ratios and the superimposed texts on the 
    grids to represent the p-values of overlaps.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\S4method{drawHeatmap}{GeneOverlapMatrix}(object, adj.p=F, cutoff=.05, 
    ncolused=9, grid.col=c("Greens", "Blues", "Greys", 
    "Oranges", "Purples", "Reds"), note.col="red")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A GeneOverlapMatrix object.}
  \item{adj.p}{Boolean label for whether p-values should be adjusted (using the 
  Benjamin-Hochberg method) before showing.}
  \item{cutoff}{P-value cutoff to mask the insignificant comparisons.}
  \item{ncolused}{Number of colors used to represent the scale of odds ratios.}
  \item{grid.col}{Color for odds ratios.}
  \item{note.col}{Color for p-value texts.}
}
\details{
    The grids that are below the p-value cutoff or have odds ratios no larger 
    than 1 will be masked and shown as the lightest color.
%%  ~~ If necessary, more details than the description above ~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(GeneOverlap)
gom.obj <- newGOM(hESC.ChIPSeq.list, genome.size=gs.RNASeq)
drawHeatmap(gom.obj, adj.p=TRUE, cutoff=1,  # show all.
    ncolused=5, grid.col="Blues", note.col="black")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
