\name{dcEnrichment}
\alias{dcEnrichment}
\title{Function to conduct enrichment analysis given the input data and the ontology in query}
\usage{
dcEnrichment(data, domain = c("SCOP.sf"), ontology = c("GOBP", "GOMF",
"GOCC"), sizeRange = c(10, 1000), min.overlap = 3,
which_distance = NULL, test = c("HypergeoTest", "FisherTest",
"BinomialTest"), p.adjust.method = c("BH", "BY", "bonferroni", "holm",
"hochberg", "hommel"), ontology.algorithm = c("none", "pc", "elim",
"lea"),
elim.pvalue = 0.01, lea.depth = 2, verbose = T,
RData.location = "http://supfam.org/dcGOR/data")
}
\arguments{
\item{data}{an input vector. It contains id for a list of
domains, for example, sunids for SCOP domains}

\item{domain}{the domain identity. It can be one of
'SCOP.sf' for SCOP superfamilies}

\item{ontology}{the ontology supported currently. It can
be "GOBP" for Gene Ontology Biological Process, "GOMF"
for Gene Ontology Molecular Function, "GOCC" for Gene
Ontology Cellular Component. For details on the
eligibility for pairs of input domain and ontology,
please refer to the online Documentations at
\url{http://supfam.org/dcGOR/docs.html}}

\item{sizeRange}{the minimum and maximum size of members
of each term in consideration. By default, it sets to a
minimum of 10 but no more than 1000}

\item{min.overlap}{the minimum number of overlaps. Only
those terms that overlap with input data at least
min.overlap (3 domains by default) will be processed}

\item{which_distance}{which distance of terms in the
ontology is used to restrict terms in consideration. By
default, it sets to 'NULL' to consider all distances}

\item{test}{the statistic test used. It can be
"FisherTest" for using fisher's exact test,
"HypergeoTest" for using hypergeometric test, or
"BinomialTest" for using binomial test. Fisher's exact
test is to test the independence between domain group
(domains belonging to a group or not) and domain
annotation (domains annotated by a term or not), and thus
compare sampling to the left part of background (after
sampling without replacement). Hypergeometric test is to
sample at random (without replacement) from the
background containing annotated and non-annotated
domains, and thus compare sampling to background. Unlike
hypergeometric test, binomial test is to sample at random
(with replacement) from the background with the constant
probability. In terms of the ease of finding the
significance, they are in order: hypergeometric test >
binomial test > fisher's exact test. In other words, in
terms of the calculated p-value, hypergeometric test <
binomial test < fisher's exact test}

\item{p.adjust.method}{the method used to adjust
p-values. It can be one of "BH", "BY", "bonferroni",
"holm", "hochberg" and "hommel". The first two methods
"BH" (widely used) and "BY" control the false discovery
rate (FDR: the expected proportion of false discoveries
amongst the rejected hypotheses); the last four methods
"bonferroni", "holm", "hochberg" and "hommel" are
designed to give strong control of the family-wise error
rate (FWER). Notes: FDR is a less stringent condition
than FWER}

\item{ontology.algorithm}{the algorithm used to account
for the hierarchy of the ontology. It can be one of
"none", "pc", "elim" and "lea". For details, please see
'Note'}

\item{elim.pvalue}{the parameter only used when
"ontology.algorithm" is "elim". It is used to control how
to declare a signficantly enriched term (and subsequently
all domains in this term are eliminated from all its
ancestors)}

\item{lea.depth}{the parameter only used when
"ontology.algorithm" is "lea". It is used to control how
many maximum depth is uded to consider the children of a
term (and subsequently all domains in these children term
are eliminated from the use for the recalculation of the
signifance at this term)}

\item{verbose}{logical to indicate whether the messages
will be displayed in the screen. By default, it sets to
TRUE for display}

\item{RData.location}{the characters to tell the location
of built-in RData files. By default, it remotely locates
at "http://supfam.org/dnet/data" or
"https://github.com/hfang-bristol/dcGOR/data". For the
user equipped with fast internet connection, this option
can be just left as default. But it is always advisable
to download these files locally. Especially when the user
needs to run this function many times, there is no need
to ask the function to remotely download every time (also
it will unnecessarily increase the runtime). For
examples, these files (as a whole or part of them) can be
first downloaded into your current working directory, and
then set this option as: \eqn{RData.location="."}. If
RData to load is already part of package itself, this
parameter can be ignored (since this function will try to
load it via function \code{data} first)}
}
\value{
an object of class "eTerm", a list with following
components: \itemize{ \item{\code{term_info}: a matrix of
nTerm X 5 containing term information, where nTerm is the
number of terms in consideration, and the 5 columns are
"term_id" (i.e. "Term ID"), "term_name" (i.e. "Term Name"),
"namespace" (i.e. "Term Namespace"), "distance" (i.e. "Term
Distance") and "IC" (i.e. "Information Content for the term
based on annotation frequency by it")} \item{\code{anno}: a
list of terms, each storing annotated domain members.
Always, terms are identified by "term_id" and domain
members identified by their ids (e.g. sunids for SCOP
domains)} \item{\code{data}: a vector containing input data
in consideration. It is not always the same as the input
data as only those mappable are retained}
\item{\code{overlap}: a list of terms, each storing domains
overlapped between domains annotated by a term and domains
in the input data (i.e. the domains of interest). Always,
terms are identified by "term_id" and domain members
identified by their ids (e.g. sunids for SCOP domains)}
\item{\code{zscore}: a vector containing z-scores}
\item{\code{pvalue}: a vector containing p-values}
\item{\code{adjp}: a vector containing adjusted p-values.
It is the p value but after being adjusted for multiple
comparisons} \item{\code{call}: the call that produced this
result} }
}
\description{
\code{dcEnrichment} is supposed to conduct enrichment
analysis given the input data and the ontology in query. It
returns an object of class "eTerm". Enrichment analysis is
based on either Fisher's exact test or Hypergeometric test.
The test can respect the hierarchy of the ontology.
}
\note{
The interpretation of the algorithms used to account for
the hierarchy of the ontology is: \itemize{ \item{"none":
does not consider the ontology hierarchy at all.}
\item{"lea": computers the significance of a term in terms
of the significance of its children at the maximum depth
(e.g. 2). Precisely, once domains are already annotated to
any children terms with a more signficance than itself,
then all these domains are eliminated from the use for the
recalculation of the signifance at that term. The final
p-values takes the maximum of the original p-value and the
recalculated p-value.} \item{"elim": computers the
significance of a term in terms of the significance of its
all children. Precisely, once domains are already annotated
to a signficantly enriched term under the cutoff of e.g.
pvalue<1e-2, all these domains are eliminated from the
ancestors of that term).} \item{"pc": requires the
significance of a term not only using the whole domains as
background but also using domains annotated to all its
direct parents/ancestors as background. The final p-value
takes the maximum of both p-values in these two
calculations.} \item{"Notes": the order of the number of
significant terms is: "none" > "lea" > "elim" > "pc".} }
}
\examples{
\dontrun{
# 1) load SCOP.sf (as 'InfoDataFrame' object)
SCOP.sf <- dcRDataLoader('SCOP.sf')
# randomly select 20 domains
data <- sample(rowNames(SCOP.sf), 20)

# 2) perform enrichment analysis
eTerm <- dcEnrichment(data, domain="SCOP.sf", ontology="GOMF")

# 3) visualise the top significant terms in the ontology hierarchy
# load obo.GOMF (as 'igraph' object)
g <- dcRDataLoader('obo.GOMF')
# focus the top 5 enriched terms
nodes_query <- names(sort(eTerm$adjp)[1:5])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
subg <- dnet::dDAGinduce(g, nodes_query)
# color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
dnet::visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]),
node.info="both", zlim=c(0,2), node.attrs=list(color=nodes.highlight))
# color-code terms according to the z-scores
dnet::visDAG(g=subg, data=eTerm$zscore[V(subg)$name], node.info="both",
colormap="darkblue-white-darkorange",
node.attrs=list(color=nodes.highlight))
}
}
\seealso{
\code{\link{dcRDataLoader}}, \code{\link{dcDAGannotate}}
}

