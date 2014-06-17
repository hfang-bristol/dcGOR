\name{SCOP.sf2GOBP}
\alias{SCOP.sf2GOBP}

\title{Annotations of SCOP domain superfamilies (sf) by Gene Ontology Biological Process (GOBP).}

\usage{
data(SCOP.sf2GOBP)
}

\description{
An object of class "Anno" that contains associations between Gene
Ontology Biological Process terms and SCOP domain superfamilies (sf).
This data is prepared based on the dcGO database (see
\url{http://supfam.org/SUPERFAMILY/dcGO/}).
}

\value{
an object of class \code{\link{Anno}}. It has slots for "annoData",
"termData" and "domainData":
\itemize{
\item{\code{annoData}: a sparse matrix of domains X terms}
\item{\code{termData}: variables describing ontology terms (i.e.
columns in annoData), including: "ID" (i.e. term ID), "Name" (i.e. term
Names), "Namespace" (i.e. term Namespace), and "Distance" (i.e. term
Distance to the ontology root)}
\item{\code{domainData}: variables describing domains (i.e. rows in
annoData), including: "id" for SCOP sunid, and "level" for SCOP level,
"description" for SCOP description}
}
}

\references{
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies
on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res},
41(Database issue):D536-44.
}

\seealso{
\code{\link{Anno-class}}
}

\keyword{datasets}

\examples{
# load data
data(SCOP.sf2GOBP)
SCOP.sf2GOBP
# retrieve info on ontology terms
termData(SCOP.sf2GOBP)
# retrieve info on SCOP domains
domainData(SCOP.sf2GOBP)
# retrieve the first 5 rows and columns of annotation data
x <- annoData(SCOP.sf2GOBP)[1:5,1:5]
x
# convert the above retrieval to the full matrix
as.matrix(x)
}
