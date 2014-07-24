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

# extract a list of domains that are present at Metazoa
flag_genome <- which(tData(Ancestral_domainome)$name=="Metazoa")
flag_domain <- annoData(Ancestral_domainome)[,flag_genome]!=0
domains_metazoa <- domainData(Ancestral_domainome)[flag_domain,]
domains_metazoa

# extract a list of domains that are present at human
flag_genome <- which(tData(Ancestral_domainome)$name=="Homo sapiens")
flag_domain <- annoData(Ancestral_domainome)[,flag_genome]!=0
domains_human <- domainData(Ancestral_domainome)[flag_domain,]
domains_human

# calculate the uniqueness and commonality
domains_metazoa_unique <- setdiff(rowNames(domains_metazoa), rowNames(domains_human))
domains_human_unique <- setdiff(rowNames(domains_human), rowNames(domains_metazoa))
domains_both <- intersect(rowNames(domains_metazoa), rowNames(domains_human))

# Enrichment analysis for domains unique in human
data <- domains_human_unique

## 1) GOMF enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOMF")
eoutput
### write into a local file <a href="GOMF_enrichments.txt">GOMF_enrichments.txt</a>
write(eoutput, file='GOMF_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput)
#### color-coded according to zscore
visEnrichment(eoutput, data.type='zscore')

## 2) GOBP enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOBP")
eoutput
### write into a local file <a href="GOBP_enrichments.txt">GOBP_enrichments.txt</a>
write(eoutput, file='GOBP_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput)

## 3) HPPA enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="HPPA")
eoutput
### write into a local file <a href="HPPA_enrichments.txt">HPPA_enrichments.txt</a>
write(eoutput, file='HPPA_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput, path.mode="all_paths", node.info="full_term_name")
#### color-coded according to zscore
visEnrichment(eoutput, data.type='zscore')

## 4) DO enrichment analysis, producing an object of S4 class 'Eoutput'
eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="DO")
eoutput
### write into a local file <a href="DO_enrichments.txt">DO_enrichments.txt</a>
write(eoutput, file='DO_enrichments.txt')
### view the top 5 significant terms
view(eoutput, top_num=5, sortBy="pvalue", details=TRUE)
### visualise the top 5 significant terms in GOMF DAG
#### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
visEnrichment(eoutput, path.mode="all_paths")

# Calculating pair-wise semantic similarity between domains unique in human
# 1) load onto.GOMF (as 'Onto' object)
g <- dcRDataLoader('onto.DO')
# 2) load SCOP superfamilies annotated by GOMF (as 'Anno' object)
Anno <- dcRDataLoader('SCOP.sf2DO')
# 3) prepare for ontology appended with annotation information
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",verbose=FALSE)
# 4) calculate pair-wise semantic similarity between domains unique in human
domains <- domains_human_unique
dnetwork <- dcDAGdomainSim(g=dag, domains=domains, method.domain="BM.average", method.term="Resnik", parallel=FALSE, verbose=TRUE)
dnetwork
# 5) convert it to an object of class 'igraph'
ig <- dcConverter(dnetwork, from='Dnetwork', to='igraph')
# 6) visualise the domain network
## extract edge weight (with 2-digit precision)
x <- signif(E(ig)$weight, digits=2)
## rescale into an interval [1,4] as edge width
edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
## do visualisation
ind <- match(V(ig)$name,domainNames(Anno))
vertex.label <- paste(V(ig)$name, '\n', as.character(dData(Anno)[ind,]$description), sep='')
dnet::visNet(g=ig, vertex.label=vertex.label, vertex.label.color="red", vertex.label.cex=0.7, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)


# RWR-based contact strength between terms
# 1) define sets of seeds as data: each seed with equal weight (i.e. all non-zero entries are '1')
Anno <- dcRDataLoader('SCOP.sf2GOMF')
flag <- match(domainNames(Anno), nodeNames(dnetwork))
ind <- which(!is.na(flag))
data <- as.matrix(annoData(Anno)[ind,])
## focus on those terms having 3 annotatable domains
ind <- apply(data,2,sum)>=3
data <- data[,ind]
# 2) calcualte their two contact graph
coutput <- dcRWRpipeline(data=data, g=dnetwork, permutation="degree", num.permutation=2000, adjp.cutoff=0.05, parallel=FALSE)
coutput
# 3) visualise the network containing contact between terms
cnet <- cnetwork(coutput)
cnet
## convert it to an object of class 'igraph'
ig <- dcConverter(cnet, from='Cnetwork', to='igraph')
## extract edge weight (with 2-digit precision)
x <- signif(E(ig)$weight, digits=2)
## rescale into an interval [1,4] as edge width
edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
## do visualisation
ind <- match(V(ig)$name,termNames(Anno))
vertex.label <- paste(V(ig)$name, '\n', tData(Anno)[ind,]$Name, sep='')
dnet::visNet(g=ig, vertex.label=vertex.label, vertex.label.color="red", vertex.label.cex=0.7, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
