\name{Pfam2GOBP}
\alias{Pfam2GOBP}

\title{Annotations of Pfam domains by Gene Ontology Biological Process (GOBP).}

\usage{
data(Pfam2GOBP)
}

\description{
An object of class "Anno" that contains associations between Gene
Ontology Biological Process terms and Pfam domains. This data is
prepared based on the Pfam database (see \url{http://pfam.xfam.org})
and \url{ftp://ftp.geneontology.org/pub/go/external2go/pfam2go}.
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
annoData), including: "id" for Pfam accession ID, and "level" always
equals "Pfam", "description" for Pfam description}
}
}

\references{
Finn et al. (2014) The Pfam protein families database. \emph{Nucleic
Acids Res}, 42(Database issue):D222-D230.
}

\seealso{
\code{\link{Anno-class}}
}

\keyword{datasets}

\examples{
# load data
data(Pfam2GOBP)
Pfam2GOBP
# retrieve info on ontology terms
termData(Pfam2GOBP)
# retrieve info on Pfam domains
domainData(Pfam2GOBP)
# retrieve the first 5 rows and columns of annotation data
x <- annoData(Pfam2GOBP)[1:5,1:5]
x
# convert the above retrieval to the full matrix
as.matrix(x)
}
