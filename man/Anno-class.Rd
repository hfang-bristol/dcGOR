% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ClassMethod-dcGOR.r
\docType{class}
\name{Anno-class}
\alias{Anno}
\alias{Anno-class}
\title{Definition for S4 class Anno}
\value{
Class Anno
}
\description{
\code{Anno} has 3 slots: annoData, termData and domainData
}
\section{Slots}{

\describe{
\item{\code{annoData}}{An object of S4 class \code{\link{AnnoData}},
containing data matrix with the column number equal to nrow(termData)
and the row number equal to nrow(domainData).}

\item{\code{termData}}{An object of S4 class
\code{\link{InfoDataFrame}}, describing information on columns in
annoData.}

\item{\code{domainData}}{An object of S4 class
\code{\link{InfoDataFrame}}, describing information on rows in
annoData.}
}}
\section{Creation}{

An object of this class can be created via: \code{new("Anno", annoData,
termData, domainData)}
}

\section{Methods}{

Class-specific methods:
\itemize{
\item{\code{dim()}: }{retrieve the dimension in the object}
\item{\code{annoData()}: }{retrieve the slot 'annoData' in the object}
\item{\code{termData()}: }{retrieve the slot 'termData' (as class
InfoDataFrame) in the object}
\item{\code{domainData()}: }{retrieve the slot 'domainData' (as class
InfoDataFrame) in the object}
\item{\code{tData()}: }{retrieve termData (as data.frame) in the
object}
\item{\code{dData()}: }{retrieve domainData (as data.frame) in the
object}
\item{\code{termNames()}: }{retrieve term names (ie, row names of
termData) in the object}
\item{\code{domanNames()}: }{retrieve domain names (ie, row names of
domainData) in the object}
}
Standard generic methods:
\itemize{
\item{\code{str()}: }{compact display of the content in the object}
\item{\code{show()}: }{abbreviated display of the object}
\item{\code{as(matrix, "Anno")}: }{convert a matrix to an object of
class Anno}
\item{\code{as(dgCMatrix, "Anno")}: }{convert a dgCMatrix (a sparse
matrix) to an object of class Anno}
\item{\code{[i,j]}: }{get the subset of the same class}
}
}

\section{Access}{

Ways to access information on this class:
\itemize{
\item{\code{showClass("Anno")}: }{show the class definition}
\item{\code{showMethods(classes="Anno")}: }{show the method definition
upon this class}
\item{\code{getSlots("Anno")}: }{get the name and class of each slot in
this class}
\item{\code{slotNames("Anno")}: }{get the name of each slot in this
class}
\item{\code{selectMethod(f, signature="Anno")}: }{retrieve the
definition code for the method 'f' defined in this class}
}
}
\examples{
# create an object of class Anno, only given a matrix
annoData <- matrix(runif(50),nrow=10,ncol=5)
as(annoData, "Anno")

# create an object of class Anno, given a matrix plus information on its columns/rows
# 1) create termData: an object of class InfoDataFrame
data <- data.frame(x=1:5, y=I(LETTERS[1:5]), row.names=paste("Term",
1:5, sep="_"))
termData <- new("InfoDataFrame", data=data)
termData
# 2) create domainData: an object of class InfoDataFrame
data <- data.frame(x=1:10, y=I(LETTERS[1:10]),
row.names=paste("Domain", 1:10, sep="_"))
domainData <- new("InfoDataFrame", data=data)
domainData
# 3) create an object of class Anno
# VERY IMPORTANT: make sure having consistent names between annoData and domainData (and termData)
annoData <- matrix(runif(50),nrow=10,ncol=5)
rownames(annoData) <- rowNames(domainData)
colnames(annoData) <- rowNames(termData)
x <- new("Anno", annoData=annoData, domainData=domainData,
termData=termData)
x
# 4) look at various methods defined on class Anno
dim(x)
annoData(x)
termData(x)
tData(x)
domainData(x)
dData(x)
termNames(x)
domainNames(x)
# 5) get the subset
x[1:3,1:2]
}
\seealso{
\code{\link{Anno-method}}
}
\keyword{S4}
\keyword{classes}

