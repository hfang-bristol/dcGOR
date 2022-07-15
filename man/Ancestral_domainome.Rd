\name{Ancestral_domainome}
\alias{Ancestral_domainome}

\title{Ancestral superfamily domain repertoires in Eukaryotes}

\usage{
data(Ancestral_domainome)
}

\description{
An object of class "Anno" that contains information about domain
repertoires (a complete set of domains: domain-ome) in Eukaryotes
(including extant and ancestral genomes). This data is prepared based
on 1) SUPERFAMILY database which provides domain and architecture
assignments to all completely sequenced genomes including eukaryotic
genomes; 2) ancestral domain architecture repertoires inferred by
applying Dollo parsimony to eukaryotic part of species tree of life
(sTOL), from which ancestral superfamily domain and architecture
repertoires at all branching points in eukaryotic evolution are
inferred. This allows us to list ancestral domain and architecture
repertoires that were present at these points. Based on the
observed/inferred domain and architecture repertoires, we also define
genome-specific plasticity potential for an individual domain as how
many different architectures (or architecture diversity) it can be
formed in an extant/ancestral genome. As a result, for each genome,
domain repertoires (domainome) are represented as a profile of states
on domains, where non-zero entry indicates a domain for which how many
different architectures have occurred in the genome.
}

\value{
an object of class \code{\link{Anno}}. It has slots for "annoData",
"termData" and "domainData":
\itemize{
\item{\code{annoData}: a sparse matrix of 2019 domains X 875
terms/genomes (including 438 extant genomes and 437 ancestral genomes),
with each entry telling how many different architectures a domain has
in a genome. Note: zero entry also means that this domain is absent in
the genome}
\item{\code{termData}: variables describing terms/genomes (i.e. columns
in annoData), including extant/ancestral genome information: "left_id"
(unique and used as internal id), "right_id" (used in combination with
"left_id" to define the post-ordered binary tree structure), "taxon_id"
(NCBI taxonomy id, if matched), "genome" (2-letter genome identifiers
used in SUPERFAMILY, if being extant), "name" (NCBI taxonomy name, if
matched), "rank" (NCBI taxonomy rank, if matched), "branchlength"
(branch length in relevance to the parent), and "common_name" (NCBI
taxonomy common name, if matched and existed)}
\item{\code{domainData}: variables describing domains (i.e. rows in
annoData), including information about domains: "sunid" for SCOP id,
"level" for SCOP level, "classification" for SCOP classification,
"description" for SCOP description}
}
}

\references{
Fang et al. (2013) A daily-updated tree of (sequenced) life as a
reference for genome research. \emph{Scientific reports}, 3:2015.
\cr
Morais et al. (2011) SUPERFAMILY 1.75 including a domain-centric gene
ontology method. \emph{Nucleic Acids Res}, 39(Database issue):D427-34.
\cr
Andreeva et al. (2008) Data growth and its impact on the SCOP database:
new developments. \emph{Nucleic Acids Res}, 36(Database issue):D419-425
}

\keyword{datasets}
\examples{
data(Ancestral_domainome)
Ancestral_domainome
# retrieve info on terms/genomes
termData(Ancestral_domainome)
# retrieve info on SCOP domains
domainData(Ancestral_domainome)
# retrieve the first 5 rows and columns of data
x <- annoData(Ancestral_domainome)[1:5,1:5]
x
# convert the above retrieval to the full matrix
as.matrix(x)
}
