\name{onto.UP}
\alias{onto.UP}

\title{UniProtKB UniPathway (UP).}

\usage{
data(onto.UP)
}

\description{
An R object that contains information on UniProtKB UniPathway terms.
These terms are organised as a direct acyclic graph (DAG), which is
further stored as an object of the class 'igraph' (see
\url{http://igraph.org/r/doc/aaa-igraph-package.html}). UniProtKB
UniPathway (UP) is a fully manually curated resource for the
representation and annotation of metabolic pathways, being used as a
controlled vocabulary for pathway annotation in UniProtKB.
}

\value{
an object of S4 class \link{Onto}. It has slots for "nodeInfo" and
"adjMatrix"
\itemize{
\item{\code{nodeInfo}: an object of S4 class \link{InfoDataFrame},
describing information on nodes/terms including: "term_id" (i.e. Term
ID), "term_name" (i.e. Term Name), "term_namespace" (i.e. Term
Namespace), and "term_distance" (i.e. Term Distance: the distance to
the root; always 0 for the root itself)}
\item{\code{adjMatrix}: an object of S4 class \link{AdjData},
containing adjacency data matrix, with rows for parent (arrow-outbound)
and columns for children (arrow-inbound)}
}
}

\references{
Morgat et al. (2006) UniPathway: a resource for the exploration and
annotation of metabolic pathways. \emph{Nucleic Acids Res}, 40:D761-9.
\cr
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\keyword{datasets}

\examples{
data(onto.UP)
onto.UP
}
