\name{onto.DO}
\alias{onto.DO}

\title{Disease Ontology (DO).}

\usage{
data(onto.DO)
}

\description{
An R object that contains information on Gene Ontology Biological
Process terms. These terms are organised as a direct acyclic graph
(DAG), which is further stored as an object of the class 'igraph' (see
\url{http://igraph.org/r/doc/aaa-igraph-package.html}). This data is
prepared based on
http://sourceforge.net/p/diseaseontology/code/HEAD/tree/trunk/HumanDO.obo.
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
Schriml et al. (2012) Disease Ontology: a backbone for disease semantic
integration. \emph{Nucleic Acids Res}, 40:D940-946.
\cr
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\keyword{datasets}

\examples{
data(onto.DO)
onto.DO
}
