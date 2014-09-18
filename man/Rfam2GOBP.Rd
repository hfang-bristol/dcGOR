\name{Rfam2GOBP}
\alias{Rfam2GOBP}

\title{Annotations of Rfam RNA families by Gene Ontology Biological Process (GOBP).}

\usage{
data(Rfam2GOBP)
}

\description{
An object of class "Anno" that contains associations between Gene
Ontology Biological Process terms and Rfam RNA families. This data is
prepared based on the Rfam database (see \url{http://rfam.xfam.org})
and \url{http://geneontology.org/external2go/rfam2go}.
}

\value{
an object of class \code{\link{Anno}}. It has slots for "annoData",
"termData" and "domainData":
\itemize{
\item{\code{annoData}: a sparse matrix of RNAs X terms}
\item{\code{termData}: variables describing ontology terms (i.e.
columns in annoData), including: "ID" (i.e. term ID), "Name" (i.e. term
Names), "Namespace" (i.e. term Namespace), and "Distance" (i.e. term
Distance to the ontology root)}
\item{\code{domainData}: variables describing RNAs (i.e. rows in
annoData), including: "id" for Rfam accession ID, and "level" always
equals "Rfam", "description" for Rfam description}
}
}

\references{
Gardner et al. (2011) Rfam: Wikipedia, clans and the "decimal" release.
\emph{Nucleic Acids Res}, 39(Database issue):D141-D145.
}

\seealso{
\code{\link{Anno-class}}
}

\keyword{datasets}

\examples{
# load data
data(Rfam2GOBP)
Rfam2GOBP
# retrieve info on ontology terms
termData(Rfam2GOBP)
# retrieve info on Rfam RNAs
domainData(Rfam2GOBP)
# retrieve the first 5 rows and columns of annotation data
x <- annoData(Rfam2GOBP)[1:5,1:5]
x
# convert the above retrieval to the full matrix
as.matrix(x)
}
