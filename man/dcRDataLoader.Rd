\name{dcRDataLoader}
\alias{dcRDataLoader}
\title{Function to load dcGOR built-in RData}
\usage{
dcRDataLoader(RData = c(NA, "SCOP.sf", "SCOP.fa", "obo.GOBP",
"obo.GOMF",
"obo.GOCC", "obo.DO", "obo.HPPA", "obo.HPMI", "obo.HPON", "obo.MP",
"obo.EC",
"obo.KW", "obo.UP", "SCOP.sf2GOBP", "SCOP.sf2GOMF", "SCOP.sf2GOCC",
"SCOP.sf2DO", "SCOP.sf2HPPA", "SCOP.sf2HPMI", "SCOP.sf2HPON",
"SCOP.sf2MP",
"SCOP.sf2EC", "SCOP.sf2KW", "SCOP.sf2UP", "SCOP.fa2GOBP",
"SCOP.fa2GOMF",
"SCOP.fa2GOCC", "SCOP.fa2DO", "SCOP.fa2HPPA", "SCOP.fa2HPMI",
"SCOP.fa2HPON",
"SCOP.fa2MP", "SCOP.fa2EC", "SCOP.fa2KW", "SCOP.fa2UP",
"Ancestral_domainome",
"eTOL"), domain = c(NA, "SCOP.sf", "SCOP.fa"), ontology = c(NA,
"GOBP", "GOMF", "GOCC", "DO", "HPPA", "HPMI", "HPON", "MP", "EC", "KW",
"UP"),
verbose = T, RData.location = "http://supfam.org/dcGOR/data")
}
\arguments{
\item{RData}{which built-in RData to load. If NOT NA,
this RData will be always loaded. It can be: domains
(including 'SCOP.sf', 'SCOP.fa'), ontologies (including
'obo.GOBP', 'obo.GOMF', 'obo.GOCC', 'obo.DO', 'obo.HPPA',
'obo.HPMI', 'obo.HPON', 'obo.MP', 'obo.EC', 'obo.KW',
'obo.UP'), annotations (including 'SCOP.sf2GOBP',
'SCOP.sf2GOMF', 'SCOP.sf2GOCC', 'SCOP.sf2DO',
'SCOP.sf2HPPA', 'SCOP.sf2HPMI', 'SCOP.sf2HPON',
'SCOP.sf2MP', 'SCOP.sf2EC', 'SCOP.sf2KW', 'SCOP.sf2UP',
'SCOP.fa2GOBP', 'SCOP.fa2GOMF', 'SCOP.fa2GOCC',
'SCOP.fa2DO', 'SCOP.fa2HPPA', 'SCOP.fa2HPMI',
'SCOP.fa2HPON', 'SCOP.fa2MP', 'SCOP.fa2EC', 'SCOP.fa2KW',
'SCOP.fa2UP'), and domainome in eukaryotic genomes
(including 'Ancestral_domainome', 'eTOL'). On the
meanings, please refer to the Documentations}

\item{domain}{domain part of annotation RData to load.
When RData is NA and this plus next are NOT NA, then this
plus next one are used to specify which annotation RData
to load. In addition to NA, it can also be: 'SCOP.sf',
'SCOP.fa'}

\item{ontology}{ontology part of annotation RData to
load. This only works together with the previous 'domain'
parameter. In addition to NA, it can also be: 'GOBP',
'GOMF', 'GOCC', 'DO', 'HPPA', 'HPMI', 'HPON', 'MP', 'EC',
'KW', 'UP'}

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
any use-specified variable that is given on the right side
of the assigement sign '<-', which contains the loaded
RData.
}
\description{
\code{dcRDataLoader} is supposed to load RData that are
used by package dcGOR.
}
\note{
If there are no use-specified variable that is given on the
right side of the assigement sign '<-', then no RData will
be loaded onto the working environment.
}
\examples{
# Always, load from specified RData directly
SCOP.sf <- dcRDataLoader(RData='SCOP.sf')
obo.GOMF <- dcRDataLoader(RData='obo.GOMF')
# But for annotaion data, there are two ways to do so:
# 1) in a direct way
SCOP.sf2GOMF <- dcRDataLoader(RData='SCOP.sf2GOMF')
# 2) in an indirect way: specify both domain and ontology
SCOP.sf2GOMF <- dcRDataLoader(domain='SCOP.sf', ontology='GOMF')
}
\seealso{
\code{\link{dcEnrichment}}
}

