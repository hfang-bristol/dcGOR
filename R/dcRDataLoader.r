#' Function to load dcGOR built-in RData
#'
#' \code{dcRDataLoader} is supposed to load RData that are used by package dcGOR.
#'
#' @param RData which built-in RData to load. If NOT NA, this RData will be always loaded. It can be: domains/RNAs (including 'SCOP.sf', 'SCOP.fa', 'Pfam', 'InterPro', 'Rfam'), ontologies (including 'onto.GOBP', 'onto.GOMF', 'onto.GOCC', 'onto.DO', 'onto.HPPA', 'onto.HPMI', 'onto.HPON', 'onto.MP', 'onto.EC', 'onto.KW', 'onto.UP'), annotations (including 'SCOP.sf2GOBP', 'SCOP.sf2GOMF', 'SCOP.sf2GOCC', 'SCOP.sf2DO', 'SCOP.sf2HPPA', 'SCOP.sf2HPMI', 'SCOP.sf2HPON', 'SCOP.sf2MP', 'SCOP.sf2EC', 'SCOP.sf2KW', 'SCOP.sf2UP', 'SCOP.fa2GOBP', 'SCOP.fa2GOMF', 'SCOP.fa2GOCC', 'SCOP.fa2DO', 'SCOP.fa2HPPA', 'SCOP.fa2HPMI', 'SCOP.fa2HPON', 'SCOP.fa2MP', 'SCOP.fa2EC', 'SCOP.fa2KW', 'SCOP.fa2UP', 'Pfam2GOBP', 'Pfam2GOMF', 'Pfam2GOCC', 'InterPro2GOBP', 'InterPro2GOMF', 'InterPro2GOCC', 'Rfam2GOBP', 'Rfam2GOMF', 'Rfam2GOCC'), domainome in eukaryotic genomes (including 'Ancestral_domainome', 'eTOL'), and databases used for predictiing ontology terms from input protein domain contents. On the meanings, please refer to the Documentations
#' @param domain domain part of annotation RData to load. When RData is NA and this plus next are NOT NA, then this plus next one are used to specify which annotation RData to load. In addition to NA, it can also be: 'SCOP.sf', 'SCOP.fa', 'Pfam' and 'InterPro'
#' @param ontology ontology part of annotation RData to load. This only works together with the previous 'domain' parameter. In addition to NA, it can also be: 'GOBP', 'GOMF', 'GOCC', 'DO', 'HPPA', 'HPMI', 'HPON', 'MP', 'EC', 'KW', 'UP'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.location the characters to tell the location of built-in RData files. By default, it remotely locates at "http://supfam.org/dcGOR/data" or "http://dcgor.r-forge.r-project.org/data". For the user equipped with fast internet connection, this option can be just left as default. But it is always advisable to download these files locally. Especially when the user needs to run this function many times, there is no need to ask the function to remotely download every time (also it will unnecessarily increase the runtime). For examples, these files (as a whole or part of them) can be first downloaded into your current working directory, and then set this option as: \eqn{RData.location="."}. If RData to load is already part of package itself, this parameter can be ignored (since this function will try to load it via function \code{data} first)
#' @return 
#' any use-specified variable that is given on the right side of the assigement sign '<-', which contains the loaded RData.
#' @note If there are no use-specified variable that is given on the right side of the assigement sign '<-', then no RData will be loaded onto the working environment. 
#' @export
#' @seealso \code{\link{dcEnrichment}}
#' @include dcRDataLoader.r
#' @examples
#' # Always, load from specified RData directly
#' SCOP.sf <- dcRDataLoader(RData='SCOP.sf')
#' Pfam <- dcRDataLoader(RData='Pfam')
#' InterPro <- dcRDataLoader(RData='InterPro')
#' Rfam <- dcRDataLoader(RData='Rfam')
#' onto.GOMF <- dcRDataLoader(RData='onto.GOMF')
#' # But for annotaion data, there are two ways to do so:
#' # 1) in a direct way
#' SCOP.sf2GOMF <- dcRDataLoader(RData='SCOP.sf2GOMF')
#' # 2) in an indirect way: specify both domain and ontology
#' SCOP.sf2GOMF <- dcRDataLoader(domain='SCOP.sf', ontology='GOMF')

dcRDataLoader <- function(RData=c(NA,'SCOP.sf','SCOP.fa','Pfam','InterPro','Rfam','onto.GOBP','onto.GOMF','onto.GOCC','onto.DO','onto.HPPA','onto.HPMI','onto.HPON','onto.MP','onto.EC','onto.KW','onto.UP','SCOP.sf2GOBP','SCOP.sf2GOMF','SCOP.sf2GOCC','SCOP.sf2DO','SCOP.sf2HPPA','SCOP.sf2HPMI','SCOP.sf2HPON','SCOP.sf2MP','SCOP.sf2EC','SCOP.sf2KW','SCOP.sf2UP','SCOP.fa2GOBP','SCOP.fa2GOMF','SCOP.fa2GOCC','SCOP.fa2DO','SCOP.fa2HPPA','SCOP.fa2HPMI','SCOP.fa2HPON','SCOP.fa2MP','SCOP.fa2EC','SCOP.fa2KW','SCOP.fa2UP','Pfam2GOBP','Pfam2GOMF','Pfam2GOCC','InterPro2GOBP','InterPro2GOMF','InterPro2GOCC','Rfam2GOBP','Rfam2GOMF','Rfam2GOCC','Ancestral_domainome','eTOL','Feature2GOBP.sf','Feature2GOMF.sf','Feature2GOCC.sf','Feature2HPPA.sf','Feature2GOBP.pfam','Feature2GOMF.pfam','Feature2GOCC.pfam','Feature2HPPA.pfam','Feature2GOBP.interpro','Feature2GOMF.interpro','Feature2GOCC.interpro','Feature2HPPA.interpro'), domain=c(NA,'SCOP.sf','SCOP.fa','Pfam','InterPro','Rfam'), ontology=c(NA,'GOBP','GOMF','GOCC','DO','HPPA','HPMI','HPON','MP','EC','KW','UP'), verbose=T, RData.location="http://dcgor.r-forge.r-project.org/data")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NA means to take the first one
    RData <- match.arg(RData)
    domain <- match.arg(domain)
    ontology <- match.arg(ontology)
    
    ###############################
    if(is.na(RData)){
        if(is.na(domain) | is.na(ontology)){
            stop(sprintf("Please make sure that either 'RData' is NOT NA or both 'domain' and 'ontology' are NOT NA.\n"))
        }else{
            
            ###############################
            ## get RData for annotations from input 'domain' and 'ontology'
            ###############################
            
            ## check the eligibility for pairs of input domain and ontology
            all.ontologies <- c('GOBP','GOMF','GOCC','DO','HPPA','HPMI','HPON','MP','EC','KW','UP')
            possible.ontologies <- switch(domain,
                               SCOP.sf = all.ontologies[c(1:11)],
                               SCOP.fa = all.ontologies[c(1:11)],
                               Pfam = all.ontologies[c(1:3)],
                               InterPro = all.ontologies[c(1:3)],
                               Rfam = all.ontologies[c(1:3)],
                               )
            if(!(ontology %in% possible.ontologies)){
                stop(sprintf("The input pair of domain (%s) and ontology (%s) are not supported.\nThe supported ontologies in domain (%s): %s.\n", domain, ontology, domain, paste(possible.ontologies,collapse=", ")))
            }else{
                RData <- paste(domain, ontology, sep='2')
            }
        }
    }
    
    ###############################

    ## make sure there is no "/" at the end
    path_host <- gsub("/$", "", RData.location)
    if(path_host=="" || length(path_host)==0 || is.na(path_host)){
        path_host <- "http://supfam.org/dcGOR/data"
    }

    ## load 
    load_remote <- paste(path_host, "/", RData, ".RData", sep="")
    load_local1 <- file.path(path_host, paste("data/", RData, ".RData", sep=""))
    load_local2 <- file.path(path_host, paste(RData, ".RData", sep=""))
    load_package <- RData

    ## first, load data from the package itself
    if(length(suppressWarnings(tryCatch(eval(parse(text=paste("data(",load_package,", package='dcGOR')",sep=""))), error=function(e) e, warning=function(w) w)))==2){
        ## second, load local R files
        RData_local <- c(load_local1, load_local2)
        load_flag <- sapply(RData_local, function(x){
            if(.Platform$OS.type=="windows") x <- gsub("/", "\\\\", x)
            ifelse(file.exists(x), TRUE, FALSE)
        })
        ## otherwise, load remote R files
        if(sum(load_flag)==0){
        
            con <- url(load_remote)
            if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                load_remote <- paste("http://supfam.org/dcGOR/data/", RData, ".RData", sep="")
            
                con <- url(load_remote)
                if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                    load_remote <- paste("http://dcgor.r-forge.r-project.org/data/", RData, ".RData", sep="")
                
                    con <- url(load_remote)
                    if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                        load_remote <- paste("https://github.com/hfang-bristol/dcGOR/tree/master/data/", RData, ".RData", sep="")
                        
                        con <- url(load_remote)
                        if(class(suppressWarnings(try(load(con), T)))=="try-error"){
                            stop("Built-in Rdata files cannot be loaded. Please check your internet connection or their location in your local machine.\n")   
                        }
                    }
                }
            }
            close(con)
        
            load_RData <- load_remote
        }else{
            load_RData <- RData_local[load_flag]
            load(load_RData)
        }
    }else{
        load_RData <- sprintf("package 'dcGOR' version %s", packageVersion("dcGOR"))
    }
    
    out <- ''
    eval(parse(text=paste("out <- ", RData, sep="")))
    
    if(verbose){
        message(sprintf("'%s' (from %s) has been loaded into the working environment", RData, load_RData), appendLF=T)
    }

    invisible(out)
}