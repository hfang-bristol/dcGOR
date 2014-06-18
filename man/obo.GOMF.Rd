\name{obo.GOMF}
\alias{obo.GOMF}

\title{Gene Ontology Molecular Function (GOMF).}

\usage{
data(obo.GOMF)
}

\description{
An R object that contains information on Gene Ontology Molecular
Function terms. These terms are organised as a direct acyclic graph
(DAG), which is further stored as an object of the class 'igraph' (see
\url{http://igraph.org/r/doc/aaa-igraph-package.html}). This data is
prepared based on
\url{http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo}.
Only the edges with the relation (either 'is_a' or 'part_of') are
retained.
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
Ashburner et al. (2000) Gene ontology: tool for the unification of
biology. \emph{Nat Genet}, 25:25-29.
\cr
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\keyword{datasets}

\examples{
data(obo.GOMF)
obo.GOMF
}
