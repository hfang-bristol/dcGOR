\name{dcConverter}
\alias{dcConverter}
\title{Function to convert an object between classes 'Onto' and 'igraph'}
\usage{
dcConverter(obj, from = c("Onto", "igraph"), to = c("igraph", "Onto"),
verbose = TRUE)
}
\arguments{
\item{obj}{an object of class "Onto" or "igraph"}

\item{from}{a character specifying the class converted
from. It can be one of "Onto" and "igraph"}

\item{to}{a character specifying the class converted to.
It can be one of "igraph" and "Onto"}

\item{verbose}{logical to indicate whether the messages
will be displayed in the screen. By default, it sets to
true for display}
}
\value{
an object of class "Onto" or "igraph"
}
\description{
\code{dcConverter} is supposed to convert an object between
classes 'Onto' and 'igraph'.
}
\note{
none
}
\examples{
# 1) load onto.GOMF (as 'Onto' object)
on <- dcRDataLoader('onto.GOMF')
on

# 2) convert the object from 'Onto' to 'igraph' class
ig <- dcConverter(on, from='Onto', to='igraph')
ig

# 3) convert the object from 'igraph' to 'Onto' class
dcConverter(ig, from='igraph', to='Onto')
}
\seealso{
\code{\link{dcRDataLoader}}, \code{\link{Onto-class}}
}

