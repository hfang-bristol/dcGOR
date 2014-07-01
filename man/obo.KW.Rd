\name{obo.KW}
\alias{obo.KW}

\title{UniProtKB KeyWords (KW).}

\usage{
data(obo.KW)
}

\description{
An R object that contains information on UniProtKB KeyWords terms.
These terms are organised as a direct acyclic graph (DAG), which is
further stored as an object of the class 'igraph' (see
\url{http://igraph.org/r/doc/aaa-igraph-package.html}). UniProtKB
KeyWords (KW) controlled vocabulary provides a summary of the entry
content and is used to index UniProtKB/Swiss-Prot entries based on 10
categories (the category "Technical term" being excluded here).
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
Bairoch et al. (2005) The Universal Protein Resource (UniProt).
\emph{Nucleic Acids Res}, 33:D154-9.
\cr
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\keyword{datasets}

\examples{
data(obo.KW)
obo.KW
}
