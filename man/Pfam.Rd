\name{Pfam}
\alias{Pfam}

\title{Pfam domains (Pfam).}

\usage{
data(Pfam)
}

\description{
An object of class "InfoDataFrame" that contains information on Pfam
domains (Pfam). This data is prepared based on
\url{ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/database_files/pfamA.txt.gz}.
}

\value{
an object of class \code{\link{InfoDataFrame}}. It has slots for data
and dimLabels:
\itemize{
\item{\code{data}: a data.frame containing information about 14831
domains (in rows), with 3 columns ("id" for Pfam accession ID, and
"level" always equals "Pfam", "description" for Pfam description)}
\item{\code{dimLabels}: a character describing labels for rows and
columns in data}
}
}

\references{
Finn et al. (2014) The Pfam protein families database. \emph{Nucleic
Acids Res}, 42(Database issue):D222-D230.
}

\seealso{
\code{\link{InfoDataFrame-class}}
}

\keyword{datasets}

\examples{
# load data
data(Pfam)
Pfam
# retrieve the dimension
dim(Pfam)
# retrieve names of columns
colNames(Pfam)
# retrieve the first 5 rows of data
Data(Pfam)[1:5,]
}
