\name{Rfam}
\alias{Rfam}

\title{RNA families (Rfam).}

\usage{
data(Rfam)
}

\description{
An object of class "InfoDataFrame" that contains information on RNA
families (Rfam). This data is prepared based on
\url{ftp://anonymous@ftp.sanger.ac.uk/pub/databases/Rfam/11.0/database_files/rfam.txt.gz}.
}

\value{
an object of class \code{\link{InfoDataFrame}}. It has slots for data
and dimLabels:
\itemize{
\item{\code{data}: a data.frame containing information about 2208 RNA
families (in rows), with 3 columns ("id" for Rfam accession ID, and
"level" always equals "Rfam", "description" for Rfam description)}
\item{\code{dimLabels}: a character describing labels for rows and
columns in data}
}
}

\references{
Gardner et al. (2011) Rfam: Wikipedia, clans and the "decimal" release.
\emph{Nucleic Acids Res}, 39(Database issue):D141-D145.
}

\seealso{
\code{\link{InfoDataFrame-class}}
}

\keyword{datasets}

\examples{
# load data
data(Rfam)
Rfam
# retrieve the dimension
dim(Rfam)
# retrieve names of columns
colNames(Rfam)
# retrieve the first 5 rows of data
Data(Rfam)[1:5,]
}
