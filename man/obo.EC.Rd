\name{obo.EC}
\alias{obo.EC}

\title{Enzyme Commission (EC).}

\usage{
data(obo.EC)
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
an object of class "igraph". As a direct graph, it has attributes to
vertices/nodes:
\itemize{
\item{\code{vertex attributes}: "name" (i.e. Term ID), "term_id" (i.e.
Term ID), "term_name" (i.e. Term Name), "term_namespace" (i.e. Term
Namespace), and "term_distance" (i.e. Term Distance: the distance to
the root; always 0 for the root itself)}
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
data(obo.EC)
obo.EC
}
