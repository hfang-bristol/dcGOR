\name{InterPro}
\alias{InterPro}

\title{InterPro domains (InterPro).}

\usage{
data(InterPro)
}

\description{
An object of class "InfoDataFrame" that contains information on
InterPro domains (InterPro). This data is prepared based on
\url{ftp://anonymous@ftp.ebi.ac.uk/pub/databases/interpro/Current/entry.list}.
}

\value{
an object of class \code{\link{InfoDataFrame}}. It has slots for data
and dimLabels:
\itemize{
\item{\code{data}: a data.frame containing information about 11638
annotatable domains (in rows), with 3 columns ("id" for InterPro ID,
and "level" always equals "InterPro", "description" for InterPro
description)}
\item{\code{dimLabels}: a character describing labels for rows and
columns in data}
}
}

\references{
Hunter et al. (2012) InterPro in 2011: new developments in the family
and domain prediction database. \emph{Nucleic Acids Res}, 40(Database
issue):D306-12.
}

\seealso{
\code{\link{InfoDataFrame-class}}
}

\keyword{datasets}

\examples{
# load data
data(InterPro)
InterPro
# retrieve the dimension
dim(InterPro)
# retrieve names of columns
colNames(InterPro)
# retrieve the first 5 rows of data
Data(InterPro)[1:5,]
}
