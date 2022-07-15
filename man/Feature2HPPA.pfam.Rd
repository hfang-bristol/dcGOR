\name{Feature2HPPA.pfam}
\alias{Feature2HPPA.pfam}

\title{Annotations and scores between Pfam domains/combinations and Human Phenotype Phenotypic Abnormality (HPPA).}

\usage{
Feature2HPPA.pfam <- dcRDataLoader(RData='Feature2HPPA.pfam')
}

\description{
An object of class "HIS" that contains annotations and scores between Human Phenotype Phenotypic Abnormality terms and Pfam domains/combinations. This data is prepared based on the dcGO database (see \url{http://supfam.org/SUPERFAMILY/dcGO/}).
}

\value{   
    an object of class \code{HIS}. It has three components for "hscore", "ic" and "slim": 
    \itemize{
        \item{\code{hscore}: a list of domain-related features, each with a term-named vector containing hypergeometric score (hscore)}
        \item{\code{ic}: a term-named vector containing information content (IC). Terms are ordered first by IC and then by longest-path level, making sure that for terms with the same IC, parental terms always come first}
        \item{\code{slim}: a list of four slims, each with a term-named vector containing information content (IC). Slim '1' for very general terms, '2' for general terms, '3' for specific terms, '4' for very specific terms}
    }   
}

\references{
Fang H and Gough J. (2013) dcGO: database of domain-centric ontologies on functions, phenotypes, diseases and more. \emph{Nucleic Acids Res}, 41(Database issue):D536-44.
\cr
Fang H and Gough J. (2013) A domain-centric solution to functional genomics via dcGO Predictor. \emph{BMC Bioinformatics}, 14(Suppl 3):S9.
}

\keyword{datasets}

\examples{
# load data
Feature2HPPA.pfam <- dcRDataLoader(RData='Feature2HPPA.pfam')
names(Feature2HPPA.pfam)
}
