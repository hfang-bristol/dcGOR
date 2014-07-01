\name{SCOP.fa2UP}
\alias{SCOP.fa2UP}

\title{Annotations of SCOP domain families (fa) by UniProtKB UniPathway (UP).}

\usage{
data(SCOP.fa2UP)
}

\description{
An object of class "Anno" that contains associations between UniProtKB
UniPathway terms and SCOP domain families (fa). This data is prepared
based on the dcGO database (see
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
data(SCOP.fa2UP)
SCOP.fa2UP
# retrieve info on ontology terms
termData(SCOP.fa2UP)
# retrieve info on SCOP domains
domainData(SCOP.fa2UP)
# retrieve the first 5 rows and columns of annotation data
x <- annoData(SCOP.fa2UP)[1:5,1:5]
x
# convert the above retrieval to the full matrix
as.matrix(x)
}
