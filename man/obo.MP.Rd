\name{obo.MP}
\alias{obo.MP}

\title{Mammalian Phenotype (MP).}

\usage{
data(obo.MP)
}

\description{
An R object that contains information on Mammalian Phenotype terms.
These terms are organised as a direct acyclic graph (DAG), which is
further stored as an object of the class 'igraph' (see
\url{http://igraph.org/r/doc/aaa-igraph-package.html}). This data is
prepared based on
\url{http://sourceforge.net/p/diseaseontology/code/HEAD/tree/trunk/HumanMP.obo}.
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
Smith et al. (2009) The Mammalian Phenotype Ontology: enabling robust
annotation and comparative analysis. \emph{Wiley Interdiscip Rev Syst
Biol Med}, 1:390-399.
\cr
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\keyword{datasets}

\examples{
data(obo.MP)
obo.MP
}
