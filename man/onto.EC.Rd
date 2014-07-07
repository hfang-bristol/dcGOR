\name{onto.EC}
\alias{onto.EC}

\title{Enzyme Commission (EC).}

\usage{
data(onto.EC)
}

\description{
An R object that contains information on Enzyme Commission terms. These
terms are organised as a direct acyclic graph (DAG), which is further
stored as an object of the class 'igraph' (see
\url{http://igraph.org/r/doc/aaa-igraph-package.html}). Enzyme
Commission (EC) is a resource focused on enzyme nomenclature, which is
a system of naming enzymes (protein catalysts) with Cross-references to
UniProt sequences. It uses four-digit EC numbers to define the reaction
catalysed. The first three digits are to define the reaction catalysed
and the fourth for a unique identifier (serial number).
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
Fleischmann et al. (2004) IntEnz, the integrated relational enzyme
database. \emph{Nucleic Acids Res}, 32:D434-7.
\cr
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\keyword{datasets}

\examples{
data(onto.EC)
onto.EC
}
