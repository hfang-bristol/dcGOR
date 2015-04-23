\name{onto.GOCC}
\alias{onto.GOCC}

\title{Gene Ontology Cellular Component (GOCC).}

\usage{
data(onto.GOCC)
}

\description{
An R object that contains information on Gene Ontology Cellular
Component terms. These terms are organised as a direct acyclic graph
(DAG), which is further stored as an object of the class 'igraph' (see
\url{http://igraph.org/r/doc/aaa-igraph-package.html}). This data is
prepared based on
\url{http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo}.
Only the edges with the relation (either 'is_a' or 'part_of') are
retained.
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
Ashburner et al. (2000) Gene ontology: tool for the unification of
biology. \emph{Nat Genet}, 25:25-29.
\cr
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\keyword{datasets}

\examples{
data(onto.GOCC)
onto.GOCC
}
