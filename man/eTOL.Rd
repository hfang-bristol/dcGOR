\name{eTOL}
\alias{eTOL}

\title{eukaryotic Tree Of Life (eTOL)}

\usage{
data(eTOL)
}

\description{
A 'phylo' object that contains information about eukaryotic part of
species tree of life (eTOL). It is a rooted binary tree. Tips represent
extant genomes. Since its reconstruction is guided under the NCBI
taxonomy, each internal node is either mapped onto a unique taxonomic
identifier or left empty (assumedly a hypothetical unknown ancestral
genome).
}

\value{
an object of class "phylo" with the following components:
\itemize{
\item{\code{Nnode}: the number of (internal) nodes}
\item{\code{tip.label}: a vector giving the names of the tips (i.e.,
"left_id" to define the post-ordered binary tree structure)}
\item{\code{node.label}: a vector giving the names of the internal
nodes (i.e., "left_id" to define the post-ordered binary tree
structure)}
\item{\code{genome_info}: a matrix of all nodes (including tips and
internal nodes) X 8, giving extant/ancestral genome information:
"left_id" (unique and used as internal id), "right_id" (used in
combination with "left_id" to define the post-ordered binary tree
structure), "taxon_id" (NCBI taxonomy id, if matched), "genome"
(2-letter genome identifiers used in SUPERFAMILY, if being extant),
"name" (NCBI taxonomy name, if matched), "rank" (NCBI taxonomy rank, if
matched), "branchlength" (branch length in relevance to the parent),
and "common_name" (NCBI taxonomy common name, if matched and existed)}
\item{\code{edge}: a two-column matrix of mode numeric where each row
represents an edge of the tree; the nodes and the tips are symbolized
with numbers; the tips are numbered 1, 2, ..., and the internal nodes
are numbered after the tips. For each row, the first column gives the
ancestor}
\item{\code{edge.length}: a numeric vector giving the lengths of the
branches given by 'edge'}
\item{\code{root.length}: a numeric value giving the length of the
branch at the root}
\item{\code{connectivity}: a matrix of internal nodes X all nodes
(including tips and internal nodes), with 1 for the presence of a
ancestor-descenant path, and 0 otherwise}
}
}

\references{
Fang et al. (2013) A daily-updated tree of (sequenced) life as a
reference for genome research. \emph{Scientific reports}, 3:2015.
}

\keyword{datasets}
\examples{
data(eTOL)
eTOL
# list all components
names(eTOL)
# extract information about the first 5 genomes
eTOL$genome_info[1:5,]
# look at the dimension of connectivity
dim(eTOL$connectivity)
\dontrun{
# visualise the connectivity matrix
Ntip <- length(eTOL$tip.label) # number of tips
Nnode <- eTOL$Nnode # number of internal nodes
data <- eTOL$connectivity
visHeatmapAdv(data, Rowv=FALSE,Colv=FALSE, zlim=c(0,1),
colormap="gray-black",
add.expr=abline(v=c(1,Ntip+1,(Ntip+Nnode+1))-0.5, col="white"),
key=FALSE, labRow=NA, labCol=NA)
}
}
