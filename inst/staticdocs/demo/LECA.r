# This is a demo for analysing ancestral superfamily domain repertoires in Eukaryotes
# 
# This ancestral domain-ome dataset (available from <a href="http://www.ncbi.nlm.nih.gov/pubmed/23778980" target="23778980">http://www.ncbi.nlm.nih.gov/pubmed/23778980</a>) is stored as an 'Anno' object (S4 class). It contains information about domain repertoires (a complete set of domains: domain-ome) in Eukaryotes (including extant and ancestral genomes):
## annoData(Ancestral_domainome): a sparse matrix of 2019 domains X 875 terms/genomes (including 438 extant genomes and 437 ancestral genomes), with each entry telling how many different architectures a domain has in a genome. Note: zero entry also means that this domain is absent in the genome
## termData(Ancestral_domainome): variables describing terms/genomes (i.e. columns in annoData), including extant/ancestral genome information: "left_id" (unique and used as internal id), "right_id" (used in combination with "left_id" to define the post-ordered binary tree structure), "taxon_id" (NCBI taxonomy id, if matched), "genome" (2-letter genome identifiers used in SUPERFAMILY, if being extant), "name" (NCBI taxonomy name, if matched), "rank" (NCBI taxonomy rank, if matched), "branchlength" (branch length in relevance to the parent), and "common_name" (NCBI taxonomy common name, if matched and existed)
## domainData: variables describing domains (i.e. rows in annoData), including information about domains: "sunid" for SCOP id, "level" for SCOP level, "classification" for SCOP classification, "description" for SCOP description
###############################################################################
library(dcGOR)

# load data as an 'Anno' object
Ancestral_domainome <- dcRDataLoader("Ancestral_domainome")
Ancestral_domainome

# extract a list of domains that are present at last eukarytoic common ancestor (LECA)
flag_genome <- which(Data(termData(Ancestral_domainome))$name=="Eukaryota")
flag_domain <- annoData(Ancestral_domainome)[,flag_genome]!=0
domains_leca <- domainData(Ancestral_domainome)[flag_domain,]
domains_leca

# extract a list of domains that are present at human
flag_genome <- which(Data(termData(Ancestral_domainome))$name=="Homo sapiens")
flag_domain <- annoData(Ancestral_domainome)[,flag_genome]!=0
domains_human <- domainData(Ancestral_domainome)[flag_domain,]
domains_human

# calculate the uniqueness and commonality
domains_leca_unique <- setdiff(rowNames(domains_leca), rowNames(domains_human))
domains_human_unique <- setdiff(rowNames(domains_human), rowNames(domains_leca))
domains_both <- intersect(rowNames(domains_leca), rowNames(domains_human))

# Enrichment analysis for domains unique in human
data <- domains_human_unique

## 1) GOMF enrichment analysis, producing an object of S4 class 'Eoutput'
eOutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOMF")
eOutput
### view the top 10 significant terms
view(eOutput, top_num=10, sortBy="pvalue", details=TRUE)
### visualise the top 10 significant terms in GOMF DAG
#### load obo.GOMF (as an 'igraph' object)
g <- dcRDataLoader('obo.GOMF')
g
#### focus on the top 10 significant terms (in terms of adjusted p-value) as nodes in query
nodes_query <- names(sort(adjp(eOutput))[1:10])
nodes_query
#### induce DAG only including nodes/terms in query
subg <- dnet::dDAGinduce(g, nodes_query)
subg
#### color-code terms according to adjusted p-values (taking the form of 10-based negative logarithm)
dnet::visDAG(g=subg, data=-1*log10(adjp(eOutput)), node.info="both")
#### also highlight (framed in black) nodes/terms in query
nodes.highlight <- rep("black", length(nodes_query))
names(nodes.highlight) <- nodes_query
dnet::visDAG(g=subg, data=-1*log10(adjp(eOutput)), node.info="both", node.attrs=list(color=nodes.highlight))
#### the same as above but color-code terms according to the z-scores
dnet::visDAG(g=subg, data=zscore(eOutput), node.info="both", colormap="darkblue-white-darkorange", node.attrs=list(color=nodes.highlight))
