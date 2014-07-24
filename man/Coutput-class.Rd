\docType{class}
\name{Coutput-class}
\alias{Coutput}
\alias{Coutput-class}
\title{Definition for S4 class Coutput}
\value{
Class Coutput
}
\description{
\code{Coutput} is an S4 class to store output by
\code{\link{dcRWRpipeline}}.
}
\section{Slots}{

\describe{
\item{\code{ratio}}{A symmetrix matrix, containing ratio}

\item{\code{zscore}}{A symmetrix matrix, containing z-scores}

\item{\code{pvalue}}{A symmetrix matrix, containing p-values}

\item{\code{adjp}}{A symmetrix matrix, containing adjusted p-values}

\item{\code{cnetwork}}{An object of S4 class \code{\link{Cnetwork}},
storing contact network.}
}}
\section{Creation}{
An object of this class can be created via:
\code{new("Coutput", ratio, zscore, pvalue, adjp,
cnetwork)}
}

\section{Methods}{
Class-specific methods: \itemize{ \item{\code{ratio()}:
}{retrieve the slot 'ratio' in the object}
\item{\code{zscore()}: }{retrieve the slot 'zscore' in
the object} \item{\code{pvalue()}: }{retrieve the slot
'pvalue' in the object} \item{\code{adjp()}: }{retrieve
the slot 'adjp' in the object} \item{\code{cnetwork()}:
}{retrieve the slot 'cnetwork' in the object}
\item{\code{write()}: }{write the object into a local
file} } Standard generic methods: \itemize{
\item{\code{str()}: }{compact display of the content in
the object} \item{\code{show()}: }{abbreviated display of
the object} }
}

\section{Access}{
Ways to access information on this class: \itemize{
\item{\code{showClass("Coutput")}: }{show the class
definition} \item{\code{showMethods(classes="Coutput")}:
}{show the method definition upon this class}
\item{\code{getSlots("Coutput")}: }{get the name and
class of each slot in this class}
\item{\code{slotNames("Coutput")}: }{get the name of each
slot in this class} \item{\code{selectMethod(f,
signature="Coutput")}: }{retrieve the definition code for
the method 'f' defined in this class} }
}
\examples{
\dontrun{
# 1) load onto.GOMF (as 'Onto' object)
g <- dcRDataLoader('onto.GOMF')

# 2) load SCOP superfamilies annotated by GOMF (as 'Anno' object)
Anno <- dcRDataLoader('SCOP.sf2GOMF')

# 3) prepare for ontology appended with annotation information
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",
verbose=TRUE)

# 4) calculate pair-wise semantic similarity between 10 randomly chosen domains
alldomains <- unique(unlist(nInfo(dag)$annotations))
domains <- sample(alldomains,10)
dnetwork <- dcDAGdomainSim(g=dag, domains=domains,
method.domain="BM.average", method.term="Resnik", parallel=FALSE,
verbose=TRUE)
dnetwork

# 5) estimate RWR dating based sample/term relationships
# define sets of seeds as data
# each seed with equal weight (i.e. all non-zero entries are '1')
data <- data.frame(aSeeds=c(1,0,1,0,1), bSeeds=c(0,0,1,0,1))
rownames(data) <- id(dnetwork)[1:5]
# calcualte their two contact graph
coutput <- dcRWRpipeline(data=data, g=dnetwork, parallel=FALSE)
coutput

# 6) write into the file 'Coutput.txt' in your local directory
write(coutput, file='Coutput.txt', saveBy="adjp")

# 7) retrieve several slots directly
ratio(coutput)
zscore(coutput)
pvalue(coutput)
adjp(coutput)
cnetwork(coutput)
}
}
\seealso{
\code{\link{Coutput-method}}
}
\keyword{S4}
\keyword{classes}

