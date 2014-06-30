library(staticdocs)
list(
    readme = "",
    
    index = list(
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
                "obo.GOCC"
            )
        ),
        sd_section("Annotations of SCOP domain superfamilies by ontologies",
            "These R objects are about SCOP domain superfamilies (sf) and their annotations by various ontologies, derived from the dcGO database.",
            c(
                "SCOP.sf",
                "SCOP.sf2GOBP",
                "SCOP.sf2GOMF",
                "SCOP.sf2GOCC"
            )
        ),
        sd_section("Complete domains (domain-ome) in Eukaryotic tree of life (eTOL)",
            "These databases are used for domain-centric genome analysis in Eukaryotes.",
            c(
                "Ancestral_domainome",
                "eTOL"
            )
        ),
        sd_section("Functions for analysis and visualisations",
            "These analysis and visualisation functions are used to process ontologies (and annotations), to do enrichment analysis, to calculate semantic similarity between annotated domains based on ontology term semantic similarity, and to perform random walk with restart upon domain-domain (semantic) networks. Most of analyses can be done via high-performance parallel computing.",
            c(
                'dcDAGannotate',
                'dcRDataLoader',
                'dcEnrichment',
                'dcDAGdomainSim',
                'dcRWRpipeline'
            )
        )

    )

)
