\docType{methods}
\name{Eoutput-method}
\alias{Eoutput-method}
\alias{adjp}
\alias{adjp,Eoutput-method}
\alias{overlap}
\alias{overlap,Eoutput-method}
\alias{pvalue}
\alias{pvalue,Eoutput-method}
\alias{show,Eoutput-method}
\alias{view}
\alias{view,Eoutput-method}
\alias{zscore}
\alias{zscore,Eoutput-method}
\title{Methods defined for S4 class Eoutput}
\usage{
overlap(x)

\S4method{overlap}{Eoutput}(x)

\S4method{zscore}{Eoutput}(x)

\S4method{pvalue}{Eoutput}(x)

\S4method{adjp}{Eoutput}(x)

\S4method{view}{Eoutput}(x, top_num = 5, sortBy = c("pvalue", "adjp",
"zscore", "nAnno", "nOverlap", "none"), decreasing = NULL, details = T)

\S4method{show}{Eoutput}(object)
}
\arguments{
\item{object}{an object of S4 class \code{Eoutput}.
Usually this is an output from
\code{\link{dcEnrichment}}}

\item{x}{an object of S4 class \code{Eoutput}. Usually
this is an output from \code{\link{dcEnrichment}}}

\item{top_num}{the maximum number (5, by default) of
terms will be viewed.}

\item{sortBy}{which statistics will be used for sorting
and viewing terms. It can be "pvalue" for p value, "adjp"
for adjusted p value, "zscore" for enrichment z-score,
"nAnno" for the number in domains annotated by a term,
"nOverlap" for the number in overlaps, and "none" for
ordering simply according to ID of terms}

\item{decreasing}{logical to indicate whether to sort in
a decreasing order. If it is null, by default it will be
true for "zscore", "nAnno" or "nOverlap"; otherwise
false}

\item{details}{logical to indicate whether the detailed
information of terms is also viewed. By default, it sets
to TRUE for the inclusion}
}
\value{
view(x) returns a data frame with following components:
\itemize{ \item{\code{term_id}: term ID}
\item{\code{nAnno}: number in domains annotated by a term}
\item{\code{nOverlap}: number in overlaps}
\item{\code{zscore}: enrichment z-score}
\item{\code{pvalue}: p value} \item{\code{adjp}: adjusted p
value} \item{\code{term_name}: term name; optional, it is
only appended when "details" is true}
\item{\code{term_namespace}: term namespace; optional, it
is only appended when "details" is true}
\item{\code{term_distance}: term distance; optional, it is
only appended when "details" is true} }
}
\description{
Methods defined for S4 class \code{Eoutput}.
}
\seealso{
\code{\link{Eoutput-class}}
}
\keyword{S4}
\keyword{methods}

