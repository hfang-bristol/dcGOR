\name{InterPro2GOMF}
\alias{InterPro2GOMF}

\title{Annotations of InterPro domains by Gene Ontology Molecular Function (GOMF).}

\usage{
data(InterPro2GOMF)
}

\description{
An object of class "Anno" that contains associations between Gene
Ontology Molecular Function terms and InterPro domains. This data is
prepared based on the InterPro database (see
\url{http://www.ebi.ac.uk/interpro/}) and
\url{ftp://anonymous@ftp.ebi.ac.uk/pub/databases/interpro/Current/interpro2go}.
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
annoData), including: "id" for InterPro ID, and "level" always equals
"InterPro", "description" for InterPro description}
}
}

\references{
Hunter et al. (2012) Manual GO annotation of predictive protein
signatures: the InterPro approach to GO curation. \emph{Database
(Oxford)}, 2012:bar068.
}

\seealso{
\code{\link{Anno-class}}
}

\keyword{datasets}

\examples{
# load data
data(InterPro2GOMF)
InterPro2GOMF
# retrieve info on ontology terms
termData(InterPro2GOMF)
# retrieve info on InterPro domains
domainData(InterPro2GOMF)
# retrieve the first 5 rows and columns of annotation data
x <- annoData(InterPro2GOMF)[1:5,1:5]
x
# convert the above retrieval to the full matrix
as.matrix(x)
}
