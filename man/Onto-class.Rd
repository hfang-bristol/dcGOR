\docType{class}
\name{Onto-class}
\alias{Onto}
\alias{Onto-class}
\title{Definition for S4 class Onto}
\value{
Class Onto
}
\description{
\code{Onto} has 2 slots: nodeInfo and adjMatrix
}
\section{Slots}{

\describe{
\item{\code{nodeInfo}}{An object of S4 class
\code{\link{InfoDataFrame}}, describing information on nodes/terms.}

\item{\code{adjMatrix}}{An object of S4 class \code{\link{AdjData}},
containing adjacency data matrix (for a direct graph), with rows for
parent (arrow-outbound) and columns for children (arrow-inbound)}
}}
\section{Creation}{
An object of this class can be created via:
\code{new("Onto", nodeInfo, adjMatrix)}
}

\section{Methods}{
Class-specific methods: \itemize{ \item{\code{dim()}:
}{retrieve the dimension in the object}
\item{\code{adjMatrix()}: }{retrieve the slot 'adjMatrix'
in the object} \item{\code{nodeInfo()}: }{retrieve the
slot 'nodeInfo' (as class InfoDataFrame) in the object}
\item{\code{nInfo()}: }{retrieve nodeInfo (as data.frame)
in the object} \item{\code{nodeNames()}: }{retrieve
node/term names (ie, row names of nodeInfo) in the
object} \item{\code{term_id()}: }{retrieve term id (ie,
column 'term_id' of nodeInfo) in the object, if any}
\item{\code{term_name()}: }{retrieve term id (ie, column
'term_name' of nodeInfo) in the object, if any}
\item{\code{term_namespace()}: }{retrieve term id (ie,
column 'term_namespace' of nodeInfo) in the object, if
any} \item{\code{term_distance()}: }{retrieve term id
(ie, column 'term_distance' of nodeInfo) in the object,
if any} } Standard generic methods: \itemize{
\item{\code{str()}: }{compact display of the content in
the object} \item{\code{show()}: }{abbreviated display of
the object} \item{\code{as(matrix, "Onto")}: }{convert a
matrix to an object of class Onto}
\item{\code{as(dgCMatrix, "Onto")}: }{convert a dgCMatrix
(a sparse matrix) to an object of class Onto}
\item{\code{[i]}: }{get the subset of the same class} }
}

\section{Access}{
Ways to access information on this class: \itemize{
\item{\code{showClass("Onto")}: }{show the class
definition} \item{\code{showMethods(classes="Onto")}:
}{show the method definition upon this class}
\item{\code{getSlots("Onto")}: }{get the name and class
of each slot in this class}
\item{\code{slotNames("Onto")}: }{get the name of each
slot in this class} \item{\code{selectMethod(f,
signature="Onto")}: }{retrieve the definition code for
the method 'f' defined in this class} }
}
\examples{
# create an object of class Onto, only given a matrix
adjM <- matrix(runif(25),nrow=5,ncol=5)
as(adjM, "Onto")

# create an object of class Onto, given a matrix plus information on nodes
# 1) create nodeI: an object of class InfoDataFrame
data <- data.frame(term_id=paste("Term", 1:5, sep="_"),
term_name=I(LETTERS[1:5]), term_namespace=rep("Namespace",5),
term_distance=1:5, row.names=paste("Term", 1:5, sep="_"))
nodeI <- new("InfoDataFrame", data=data)
nodeI
# 2) create an object of class Onto
# VERY IMPORTANT: make sure having consistent names between nodeInfo and adjMatrix
adjM <- matrix(runif(25),nrow=5,ncol=5)
colnames(adjM) <- rownames(adjM) <- rowNames(nodeI)
x <- new("Onto", adjMatrix=adjM, nodeInfo=nodeI)
x
# 3) look at various methods defined on class Onto
dim(x)
adjMatrix(x)
nodeInfo(x)
nInfo(x)
nodeNames(x)
term_id(x)
term_namespace(x)
term_distance(x)
# 4) get the subset
x[1:2]
}
\seealso{
\code{\link{Onto-method}}
}
\keyword{S4}
\keyword{classes}

