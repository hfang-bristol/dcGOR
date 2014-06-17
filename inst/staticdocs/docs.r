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
        sd_section("Annotations of SCOP domain superfamilies by ontologies",
            "These R objects are about SCOP domain superfamilies (sf) and their annotations by various ontologies, derived from the dcGO database.",
            c(
                "SCOP.sf",
                "SCOP.sf2GOBP",
                "SCOP.sf2GOMF",
                "SCOP.sf2GOCC"
            )
        )
    )

)
