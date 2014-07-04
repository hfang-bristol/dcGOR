\docType{class}
\name{InfoDataFrame-class}
\alias{InfoDataFrame}
\alias{InfoDataFrame-class}
\title{Definition for S4 class InfoDataFrame}
\value{
Class InfoDataFrame
}
\description{
\code{InfoDataFrame} has two slots: data and dimLabels.
}
\section{Slots}{

\describe{
\item{\code{data}}{A data.frame containing terms (rows) and measured
variables (columns).}

\item{\code{dimLabels}}{A character descripting labels for rows and
columns.}
}}
\section{Creation}{
An object of this class can be created via:
\code{new("InfoDataFrame", data, dimLabels)}
}

\section{Methods}{
Class-specific methods: \itemize{ \item{\code{dim()}:
}{retrieve the dimension in the object}
\item{\code{nrow()}: }{retrieve number of rows in the
object} \item{\code{ncol()}: }{retrieve number of columns
in the object} \item{\code{rowNames()}: }{retrieve names
of rows in the object} \item{\code{colNames()}:
}{retrieve names of columns in the object}
\item{\code{dimLabels()}: }{retrieve the slot
'dimLabels', containing labels used for display of rows
and columns in the object} \item{\code{Data()}:
}{retrieve the slot 'data' in the object} } Standard
generic methods: \itemize{ \item{\code{str()}: }{compact
display of the content in the object}
\item{\code{show()}: }{abbreviated display of the object}
\item{\code{as(data.frame, "InfoDataFrame")}: }{convert a
data.frame to an object of class InfoDataFrame}
\item{\code{[i,j]}: }{get the subset of the same class} }
}

\section{Access}{
Ways to access information on this class: \itemize{
\item{\code{showClass("InfoDataFrame")}: }{show the class
definition}
\item{\code{showMethods(classes="InfoDataFrame")}: }{show
the method definition upon this class}
\item{\code{getSlots("InfoDataFrame")}: }{get the name
and class of each slot in this class}
\item{\code{slotNames("InfoDataFrame")}: }{get the name
of each slot in this class} \item{\code{selectMethod(f,
signature="InfoDataFrame")}: }{retrieve the definition
code for the method 'f' defined in this class} }
}
\examples{
data <- data.frame(x=1:5, y=I(LETTERS[1:5]), row.names=paste("Term",
1:5, sep="_"))
dimLabels <- c("rowLabels", "colLabels")
# create an object of class InfoDataFrame
x <- new("InfoDataFrame", data=data, dimLabels=dimLabels)
x
# alternatively, using coerce methods
x <- as(data, "InfoDataFrame")
x
# look at various methods defined on class Anno
dimLabels(x)
dim(x)
nrow(x)
ncol(x)
rowNames(x)
colNames(x)
Data(x)
x[1:3,]
}
\seealso{
\code{\link{InfoDataFrame-method}}
}
\keyword{S4}
\keyword{classes}

