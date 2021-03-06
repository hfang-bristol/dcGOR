% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/dcSubtreeTips.r
\name{dcSubtreeTips}
\alias{dcSubtreeTips}
\title{Function to extract a tip-induced subtree from a phylo-formatted phylogenetic tree}
\usage{
dcSubtreeTips(phy, choose.tip.labels = NULL, subtree.type = c("clade",
"tips_only"), verbose = T)
}
\arguments{
\item{phy}{an object of class 'phylo'}

\item{choose.tip.labels}{a character specifying which tips are chosen}

\item{subtree.type}{a character specifying how to extract subtree from
given tips. It can be 'clade' or 'tips_only'. The former is first
induce clade (an internal node) from tip labels, and then to extract
the subtree under the induced clade. The latter is to directly extract
the subtree only containing given tip labels, (if necessary), after
trimming out unnecessary internal nodes}

\item{verbose}{logical to indicate whether the messages will be
displayed in the screen. By default, it sets to TRUE for display}
}
\value{
an object of class 'phylo'
}
\description{
\code{dcSubtreeTips} is supposed to extract a tip-induced subtree from
a phylo-formatted phylogenetic tree. In addition to the tree in
subject, another input is a vector containing tip labels of interest.
From valid tip lables, there are two types of subtree to extract. One
is first induce clade (an internal node) from tip labels, and then the
subtree is extracted under the induced clade. Another type is to
extract a subtree only containing given tip labels; in this situation,
some internal nodes perhaps need to further trimmed. The resulting
subtree is also represented as an object of class 'phylo'.
}
\note{
nonde
}
\examples{
# 1) with internal node labels
tree <- "(((t1:5,t2:5)i3:2,(t3:4,t4:4)i4:3)i2:2,(t5:4,t6:4)i5:6)i1;"
phy <- ape::read.tree(text=tree)
ape::plot.phylo(phy, type="p", use.edge.length=TRUE,
show.node.label=TRUE)

# 2) tip labels of interest
choose.tip.labels <- c('t1','t2','t3')
# 2a) extract subtree via an induced clade
subphy <- dcSubtreeTips(phy, choose.tip.labels, subtree.type="clade")
ape::plot.phylo(subphy, type="p", use.edge.length=TRUE,
show.node.label=TRUE)
# 2b) extract subtree containing only tips
subphy <- dcSubtreeTips(phy, choose.tip.labels,
subtree.type="tips_only")
ape::plot.phylo(subphy, type="p", use.edge.length=TRUE,
show.node.label=TRUE)
}
\seealso{
\code{\link{dcTreeConnectivity}}, \code{\link{dcSubtreeClade}}
}

