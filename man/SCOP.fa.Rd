\name{SCOP.fa}
\alias{SCOP.fa}

\title{SCOP domain families (fa).}

\usage{
data(SCOP.fa)
}

\description{
An object of class "InfoDataFrame" that contains information on SCOP
domain families (fa).
}

\value{
an object of class \code{\link{InfoDataFrame}}. It has slots for data
and dimLabels:
\itemize{
\item{\code{data}: a data.frame containing information about 2223
domains (in rows), with 3 columns ("id" for SCOP sunid, and "level" for
SCOP level, "description" for SCOP description)}
\item{\code{dimLabels}: a character describing labels for rows and
columns in data}
}
}

\references{
Morais et al. (2011) SUPERFAMILY 1.75 including a domain-centric gene
ontology method. \emph{Nucleic Acids Res}, 39(Database issue):D427-34.
\cr
Andreeva et al. (2008) Data growth and its impact on the SCOP database:
new developments. \emph{Nucleic Acids Res}, 36(Database issue):D419-425
}

\seealso{
\code{\link{InfoDataFrame-class}}
}

\keyword{datasets}

\examples{
# load data
data(SCOP.fa)
SCOP.fa
# retrieve the dimension
dim(SCOP.fa)
# retrieve names of columns
colNames(SCOP.fa)
# retrieve the first 5 rows of data
Data(SCOP.fa)[1:5,]
}
