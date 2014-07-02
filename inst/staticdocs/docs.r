library(staticdocs)
list(
    readme = "",
    
    index = list(
        sd_section("Functions for analysis and visualisations",
            "These analysis and visualisation functions are used to process ontologies (and annotations), to do enrichment analysis, to calculate semantic similarity between annotated domains based on ontology term semantic similarity, and to perform random walk with restart upon domain-domain (semantic) networks. Most of analyses can be done via high-performance parallel computing.",
            c(
                'dcDAGannotate',
                'dcRDataLoader',
                'dcEnrichment',
                'dcDAGdomainSim',
                'dcRWRpipeline'
            )
        ),
        sd_section("Definitions for S4 classes and methods",
            "These documentations are to help understand S4 classes and methods defined in the package.",
            c(
                'InfoDataFrame-class',
                'InfoDataFrame-method',
                'AnnoData-class',
                'Anno-class',
                'Anno-method'
            )
        ),
        sd_section("Ontologies mainly including open biomedical ontology (obo)",
            "These ontologies each are represented as a direct acyclic graph (DAG). DAG is stored as an object of class 'igraph'.",
            c(
                "obo.GOBP",
                "obo.GOMF",
                "obo.GOCC",
                "obo.DO",
                "obo.HPPA",
                "obo.HPMI",
                "obo.HPON",
                "obo.MP",
                "obo.EC",
                "obo.KW",
                "obo.UP"
            )
        ),
        sd_section("SCOP domain superfamilies and their annotations by ontologies",
            "These R objects are about SCOP domain superfamilies (sf) and their annotations by various ontologies, derived from the dcGO database.",
            c(
                "SCOP.sf",
                "SCOP.sf2GOBP",
                "SCOP.sf2GOMF",
                "SCOP.sf2GOCC",
                "SCOP.sf2DO",
                "SCOP.sf2HPPA",
                "SCOP.sf2HPMI",
                "SCOP.sf2HPON", 
                "SCOP.sf2MP",
                "SCOP.sf2EC",
                "SCOP.sf2KW",
                "SCOP.sf2UP"
            )
        ),
        sd_section("SCOP domain families and their annotations by ontologies",
            "These R objects are about SCOP domain families (fa) and their annotations by various ontologies, derived from the dcGO database.",
            c(
                "SCOP.fa",
                "SCOP.fa2GOBP",
                "SCOP.fa2GOMF",
                "SCOP.fa2GOCC",
                "SCOP.fa2DO",
                "SCOP.fa2HPPA",
                "SCOP.fa2HPMI",
                "SCOP.fa2HPON", 
                "SCOP.fa2MP",
                "SCOP.fa2EC",
                "SCOP.fa2KW",
                "SCOP.fa2UP"
            )
        ),
        sd_section("Complete domains (domain-ome) in Eukaryotic tree of life (eTOL)",
            "These databases are used for domain-centric genome analysis in Eukaryotes. Note: these domains are defined as SCOP domain superfamilies.",
            c(
                "Ancestral_domainome",
                "eTOL"
            )
        )

    )

)
