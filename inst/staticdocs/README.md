<a href="index.html"><IMG src="dcGOR_logo.png" height="100px" id="logo"></a>

<B><h4>An open-source <a href="http://www.r-project.org" target="R" style="font-size: 12px; color: #4169E1; text-decoration: overline; border-bottom: 1px solid #4169E1">R</a> package tailored to the needs of analysing domain-centric ontologies and annotations (<a href="http://supfam.org/SUPERFAMILY/dcGO" target="dcGO" style="font-size: 12px; color: #F87217; text-decoration: overline; border-bottom: 1px solid #F87217">dcGO</a>)</h4></B>


## Features

* `Database`: R package providing domain-centric annotations by organism-independent ontologies (eg "Gene Ontology") and organism-specific ontologies (eg "Human Phenotype" and "Mammalian Phenotype"). See <a href="docs.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">Documentations</a>
* `Infrastructure`: data structure storing ontologies (as objects of S4 class <a href="Onto-class.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">Onto</a>), annotations (as objects of S4 class <a href="Anno-class.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">Anno</a>), and enrichment outputs (as objects of S4 class <a href="Eoutput-class.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">Eoutput</a>)
* `True-path rule`: able to propagate annotations to the root. See <a href="dcDAGannotate.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">dcDAGannotate</a>
* `Enrichment analysis`: domain-based enrichment analysis and visualisation. See <a href="dcEnrichment.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">dcEnrichment</a> and <a href="visEnrichment.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">visEnrichment</a>
* `Semantic similarity`: domain-domain semantic similarity according to their annotations by an ontology. See <a href="dcDAGdomainSim.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">dcDAGdomainSim</a>
* `Random Walk with Restart`: support for walk on domain semantic similarity network.  See <a href="dcRWRpipeline.html" style="font-size: 12px; color: #000000; text-decoration: none; border-bottom: 1px solid #000000">dcRWRpipeline</a>
* `Parallel computing`: most of analyses are supported with parallel option to reduce runtime
