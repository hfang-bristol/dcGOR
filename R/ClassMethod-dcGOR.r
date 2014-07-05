require(Matrix)
################################################################################
################################################################################
#' @title Definition for S4 class InfoDataFrame
#' @description \code{InfoDataFrame} has two slots: data and dimLabels.
#' @return Class InfoDataFrame
#' @slot data A data.frame containing terms (rows) and measured variables (columns).
#' @slot dimLabels A character descripting labels for rows and columns.
#' @section Creation:
#' An object of this class can be created via: \code{new("InfoDataFrame", data, dimLabels)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{nrow()}: }{retrieve number of rows in the object}
#' \item{\code{ncol()}: }{retrieve number of columns in the object}
#' \item{\code{rowNames()}: }{retrieve names of rows in the object}
#' \item{\code{colNames()}: }{retrieve names of columns in the object}
#' \item{\code{dimLabels()}: }{retrieve the slot 'dimLabels', containing labels used for display of rows and columns in the object}
#' \item{\code{Data()}: }{retrieve the slot 'data' in the object}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' \item{\code{as(data.frame, "InfoDataFrame")}: }{convert a data.frame to an object of class InfoDataFrame}
#' \item{\code{[i,j]}: }{get the subset of the same class}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("InfoDataFrame")}: }{show the class definition}
#' \item{\code{showMethods(classes="InfoDataFrame")}: }{show the method definition upon this class}
#' \item{\code{getSlots("InfoDataFrame")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("InfoDataFrame")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="InfoDataFrame")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @import Matrix
#' @import igraph
#' @docType class
#' @keywords S4 classes
#' @name InfoDataFrame-class
#' @rdname InfoDataFrame-class
#' @seealso \code{\link{InfoDataFrame-method}}
#' @examples
#' data <- data.frame(x=1:5, y=I(LETTERS[1:5]), row.names=paste("Term", 1:5, sep="_"))
#' dimLabels <- c("rowLabels", "colLabels")
#' # create an object of class InfoDataFrame
#' x <- new("InfoDataFrame", data=data, dimLabels=dimLabels)
#' x
#' # alternatively, using coerce methods
#' x <- as(data, "InfoDataFrame")
#' x
#' # look at various methods defined on class Anno
#' dimLabels(x)
#' dim(x)
#' nrow(x)
#' ncol(x)
#' rowNames(x)
#' colNames(x)
#' Data(x)
#' x[1:3,]

#' @rdname InfoDataFrame-class
#' @aliases InfoDataFrame
#' @exportClass InfoDataFrame
setClass(
    Class="InfoDataFrame",
    representation(
        data = "data.frame",
        dimLabels = "character"
    ),
    prototype = prototype(
        data = new( "data.frame" ),
        dimLabels=c("rowLabels", "colLabels")
    ),
    validity = function(object){
        if(!is.data.frame(object@data)){
            return("data is not data.frame")
        }else{
            return(TRUE)
        }
    } 
)

########################################
#' @title Methods defined for S4 class InfoDataFrame
#' @description Methods defined for class \code{InfoDataFrame}.
#' @param x an object of class \code{InfoDataFrame}
#' @param object an object of class \code{InfoDataFrame}
#' @param i an index
#' @param j an index
#' @param ... additional parameters
#' @docType methods
#' @keywords S4 methods
#' @name InfoDataFrame-method
#' @rdname InfoDataFrame-method
#' @seealso \code{\link{InfoDataFrame-class}}

setGeneric("dimLabels", function(x) standardGeneric("dimLabels"))
#' @rdname InfoDataFrame-method
#' @aliases dimLabels
#' @export
setMethod("dimLabels", "InfoDataFrame", function(x) x@dimLabels)

#' @rdname InfoDataFrame-method
#' @aliases dim,InfoDataFrame-method
#' @export
setMethod("dim", "InfoDataFrame", function(x) dim(x@data))

#' @rdname InfoDataFrame-method
#' @aliases nrow
#' @export
setMethod("nrow", "InfoDataFrame", function(x) nrow(x@data))

#' @rdname InfoDataFrame-method
#' @aliases ncol
#' @export
setMethod("ncol", "InfoDataFrame", function(x) ncol(x@data))

setGeneric("rowNames", function(x) standardGeneric("rowNames"))
#' @rdname InfoDataFrame-method
#' @aliases rowNames
#' @export
setMethod("rowNames", "InfoDataFrame", function(x) rownames(x@data))

setGeneric("colNames", function(x) standardGeneric("colNames"))
#' @rdname InfoDataFrame-method
#' @aliases colNames
#' @export
setMethod("colNames", "InfoDataFrame", function(x) colnames(x@data))

setGeneric("Data", function(x) standardGeneric("Data"))
#' @rdname InfoDataFrame-method
#' @aliases Data
#' @export
setMethod("Data", "InfoDataFrame", function(x) x@data)

#' @rdname InfoDataFrame-method
#' @name data.frame2InfoDataFrame
setAs("data.frame", "InfoDataFrame", function(from) new("InfoDataFrame",data=from))

.wrapcat <- function(lbl, nms, total, ..., indent=2, exdent=4)
{
    lbl <- sprintf("%s:", lbl)
    txt <- paste(c(lbl,  nms), collapse=" ")
    ext <-
        if (length(nms) < total) sprintf("(%d total)", total)
        else character()
    txt <- paste(c(lbl,  nms, ext), collapse=" ")
    cat(strwrap(txt, ..., indent=indent, exdent=exdent), sep="\n")
}
.selectSomeIndex <- function(x, maxToShow=5, byrow=TRUE, ...)
{
    len <-
        if (byrow) dim(x)[[1]]
        else dim(x)[[2]]
    if (maxToShow < 3) maxToShow <- 3
    if (len > maxToShow) {
        maxToShow <- maxToShow - 1
        bot <- ceiling(maxToShow/2)
        top <- len-(maxToShow-bot-1)
        list(1:bot, "...", top:len)
    } else if (len >= 1) {
        list(1:len, NULL, NULL)
    }else{
        list(NULL, NULL, NULL)
    }
}
.showInfoDataFrame <- function(object, labels=list(0)) 
{
    lbls <- list(
        object=class(object),
        termNames=dimLabels(object)[[1]],
        varLabels="varLabels"
    )
    lbls[names(labels)] <- labels
    ## create a simplified object for extracting names
    idx <- .selectSomeIndex(Data(object), maxToShow=4)
    idy <- .selectSomeIndex(Data(object), byrow=FALSE, maxToShow=4)
    Data <- Data(object)[c(idx[[1]], idx[[3]]), c(idy[[1]], idy[[3]]), drop=FALSE]
    rnms <- rownames(Data)
    nms <- c(rnms[idx[[1]]], idx[[2]], if (!is.null(idx[[1]])) rnms[-idx[[1]]] else NULL)
    ## for terms
    .wrapcat(lbls$termNames, nms, nrow(object))
    ## for domains
    cnms <- colnames(Data)
    vars <- c(cnms[idy[[1]]], idy[[2]], cnms[-idy[[1]]])
    .wrapcat(lbls$varLabels, vars, ncol(object))
}
#' @rdname InfoDataFrame-method
#' @export
setMethod("show", "InfoDataFrame",
    function(object) {
        cat("An object of S4 class '", class(object), "'\n", sep="")
        if(sum(dim(object@data)) !=0 ){
            .showInfoDataFrame(
                object, 
                labels=list(
                    object="InfoDataFrame",
                    termNames="rowNames",
                    varLabels="colNames"
                )
            )
        }else{
            cat("but contains no data\n", sep="")
        }
    }
)

#' @rdname InfoDataFrame-method
#' @aliases [,InfoDataFrame-method
#' @export
setMethod("[", signature(x="InfoDataFrame"),
    function(x, i, j, ..., drop = FALSE) {
        if (missing(drop)){
            drop = FALSE
        }
        
        if(missing(j)) {
            labels <- x@dimLabels
            D <- x@data[i,,drop = drop]
        } else {
            labels <- x@dimLabels[j,,drop = drop]
        }
        
        if( missing( i )){
            D <- x@data[,j,drop = drop]
        }else{
            D <- x@data[i,j,drop = drop]
        }
        
        x <- new("InfoDataFrame", data=D, dimLabels=labels)
    }
)

################################################################################
################################################################################
#' @title Definition for VIRTUAL S4 class AnnoData
#' @description \code{AnnoData} is union of other classes: either matrix, or data.frame or dgCMatrix (a sparse matrix in the package Matrix). It is used as a virtual class
#' @return Class AnnoData
#' @import Matrix
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name AnnoData-class
#' @rdname AnnoData-class
#' @seealso \code{\link{Anno-class}}

#' @rdname AnnoData-class
#' @aliases AnnoData
#' @exportClass AnnoData
setClassUnion("AnnoData", c("matrix", "data.frame", "dgCMatrix"))

################################################################################
################################################################################
#' @title Definition for S4 class Anno
#' @description \code{Anno} has 3 slots: annoData, termData and domainData
#' @return Class Anno
#' @slot annoData An object of class AnnoData, containing data matrix with the column number equal to nrow(termData) and the row number equal to nrow(domainData).
#' @slot termData An object of class InfoDataFrame, describing information on columns in annoData.
#' @slot domainData An object of class InfoDataFrame, describing information on rows in annoData.
#' @section Creation:
#' An object of this class can be created via: \code{new("Anno", annoData, termData, domainData)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{annoData()}: }{retrieve the slot 'annoData' in the object}
#' \item{\code{termData()}: }{retrieve the slot 'termData' (as class InfoDataFrame) in the object}
#' \item{\code{domainData()}: }{retrieve the slot 'domainData' (as class InfoDataFrame) in the object}
#' \item{\code{tData()}: }{retrieve termData (as data.frame) in the object}
#' \item{\code{dData()}: }{retrieve domainData (as data.frame) in the object}
#' \item{\code{termNames()}: }{retrieve term names (ie, row names of termData) in the object}
#' \item{\code{domanNames()}: }{retrieve domain names (ie, row names of domainData) in the object}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' \item{\code{as(matrix, "Anno")}: }{convert a matrix (or data.frame) to an object of class Anno}
#' \item{\code{[i,j]}: }{get the subset of the same class}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Anno")}: }{show the class definition}
#' \item{\code{showMethods(classes="Anno")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Anno")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Anno")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Anno")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Anno-class
#' @seealso \code{\link{Anno-method}}
#' @examples
#' # create an object of class Anno, only given a matrix
#' annoData <- matrix(runif(50),nrow=10,ncol=5)
#' as(annoData, "Anno")
#'
#' # create an object of class Anno, given a matrix plus information on its columns/rows
#' # 1) create termData: an object of class InfoDataFrame
#' data <- data.frame(x=1:5, y=I(LETTERS[1:5]), row.names=paste("Term", 1:5, sep="_"))
#' termData <- new("InfoDataFrame", data=data)
#' termData
#' # 2) create domainData: an object of class InfoDataFrame
#' data <- data.frame(x=1:10, y=I(LETTERS[1:10]), row.names=paste("Domain", 1:10, sep="_"))
#' domainData <- new("InfoDataFrame", data=data)
#' domainData
#' # 3) create an object of class Anno
#' # VERY IMPORTANT: make sure having consistent names between annoData and domainData (and termData)
#' annoData <- matrix(runif(50),nrow=10,ncol=5)
#' rownames(annoData) <- rowNames(domainData)
#' colnames(annoData) <- rowNames(termData)
#' x <- new("Anno", annoData=annoData, domainData=domainData, termData=termData)
#' x
#' # 4) look at various methods defined on class Anno
#' dim(x)
#' annoData(x)
#' termData(x)
#' tData(x)
#' domainData(x)
#' dData(x)
#' termNames(x)
#' domainNames(x)
#' # 5) get the subset
#' x[1:3,1:2]

#' @rdname Anno-class
#' @aliases Anno
#' @exportClass Anno
setClass(
    Class="Anno",
    representation(
        annoData = "AnnoData",
        termData = "InfoDataFrame",
        domainData = "InfoDataFrame"
    ),
    prototype = prototype(
        annoData = matrix(),
        termData = new("InfoDataFrame",dimLabels=c("termNames", "termColumns")),
        domainData = new("InfoDataFrame",dimLabels=c("domainNames", "domainColumns"))
    ),
    validity = function(object){
        msg <- NULL
        # dimension for annoData
        adim <- dim(object)
        ## annoData and domainData
        if( dim(domainData(object))[1] != 0 ){
            if (adim[1] != dim(domainData(object))[1]){
                msg <- append(msg, "domain numbers differ between annoData and domainData")
            }
            if (!identical(domainNames(object), rowNames(domainData(object)))){
                msg <- append(msg, "domain names differ between annoData and domainData")
            }
        }
        ## annoData and termData
        if( dim(termData(object))[1] != 0 ){
            if (adim[2] != dim(termData(object))[1]){
                msg <- append(msg, "term numbers differ between annoData and termData")
            }
            if (!identical(termNames(object), rowNames(termData(object)))){
                msg <- append(msg, "term names differ between annoData and termData")
            }
        }
        if (is.null(msg)) TRUE else msg
    }
)

########################################
#' @title Methods defined for S4 class Anno
#' @description Methods defined for class \code{Anno}.
#' @param x an object of class \code{Anno}
#' @param object an object of class \code{Anno}
#' @param i an index
#' @param j an index
#' @param ... additional parameters
#' @docType methods
#' @keywords S4 methods
#' @name Anno-method
#' @rdname Anno-method
#' @seealso \code{\link{Anno-class}}

#' @rdname Anno-method
#' @aliases dim,Anno-method
#' @export
setMethod("dim", "Anno", function(x) dim(x@annoData))

setGeneric("annoData", function(x) standardGeneric("annoData"))
#' @rdname Anno-method
#' @aliases annoData
#' @export
setMethod("annoData", "Anno", function(x) x@annoData)

setGeneric("termData", function(x) standardGeneric("termData"))
#' @rdname Anno-method
#' @aliases termData
#' @export
setMethod("termData", "Anno", function(x) x@termData)

setGeneric("domainData", function(x) standardGeneric("domainData"))
#' @rdname Anno-method
#' @aliases domainData
#' @export
setMethod("domainData", "Anno", function(x) x@domainData)

setGeneric("tData", function(object) standardGeneric("tData"))
#' @rdname Anno-method
#' @aliases tData
#' @export
setMethod("tData", signature(object="Anno"), function(object){
    data <- Data(termData(object))
    if(sum(dim(data))==0){
        cat("No data is available\n", sep="")
    }else{
        data
    }
})

setGeneric("dData", function(object) standardGeneric("dData"))
#' @rdname Anno-method
#' @aliases dData
#' @export
setMethod("dData", signature(object="Anno"), function(object){
    data <- Data(domainData(object))
    if(sum(dim(data))==0){
        cat("No data is available\n", sep="")
    }else{
        data
    }
})

setGeneric("termNames", function(object) standardGeneric("termNames"))
#' @rdname Anno-method
#' @aliases termNames
#' @export
setMethod("termNames", signature(object="Anno"), function(object) rowNames(termData(object)))

setGeneric("domainNames", function(object) standardGeneric("domainNames"))
#' @rdname Anno-method
#' @aliases domainNames
#' @export
setMethod("domainNames", signature(object="Anno"), function(object) rowNames(domainData(object)))

#' @rdname Anno-method
#' @name matrix2Anno
setAs("matrix", "Anno", function(from) {
    ## for domainData
    rn <- rownames(from)
    if(is.null(rn)) rn <- 1:nrow(from)
    domainData <- new("InfoDataFrame", data=data.frame(names=rn))
    ## for termData    
    cn <- colnames(from)
    if(is.null(cn)) cn <- 1:ncol(from)
    termData <- new("InfoDataFrame", data=data.frame(names=cn))
    ## for Anno
    new("Anno", annoData=from, domainData=domainData, termData=termData)
})

#' @rdname Anno-method
#' @export
setMethod("show", 
    signature=signature(object="Anno"),
    function(object) {
        cat("An object of S4 class '", class(object), "'\n", sep="")
        adim <- dim(object)
        if (length(adim)>1){
            cat("annoData:", if (length(adim)>1) paste(adim[[1]], "domains,",adim[[2]], "terms") else NULL, "\n")
        }
        ## termData
        if( dim(termData(object))[1] != 0 ){
            cat("termData (", class(termData(object)), ")\n", sep="")
            .showInfoDataFrame(
                termData(object), 
                labels=list(
                    object="termData",
                    termNames="termNames",
                    varLabels="tvarLabels"
                )
            )
        }else{
            cat("termData (NULL)\n", sep="")
        }
        ## domainData
        if( dim(domainData(object))[1] != 0 ){
            cat("domainData (", class(domainData(object)), ")\n", sep="")
            .showInfoDataFrame(
                domainData(object), 
                labels=list(
                    object="domainData",
                    termNames="domainNames",
                    varLabels="dvarLabels"
                )
            )
        }else{
            cat("domainData (NULL)\n", sep="")
        }
    }
)

#' @rdname Anno-method
#' @aliases [,Anno-method
#' @export
setMethod("[", signature(x="Anno"), 
    function(x, i, j, ..., drop = FALSE) {
        if (missing(drop)){
            drop <- FALSE
        }
        if (missing(i) && missing(j)) {
            if (length(list(...))!=0){
                stop("specify domains or terms to subset")
            }
            return(x)
        }

        if (!missing(j)) {
            tD <- termData(x)[j,, ..., drop=drop]
        }else{
            tD <- termData(x)
        }
        
        if (!missing(i)) {
            dD <- domainData(x)[i,,..., drop=drop]
        }else{
            dD <- domainData(x)
        }
        
        if (missing(i) & !missing(j)){
            aD <- annoData(x)[,j]
        }else if (!missing(i) & missing(j)){
            aD <- annoData(x)[i,]
        }else if (!missing(i) & !missing(j)){
            aD <- annoData(x)[i,j]
        }else{
            aD <- annoData(x)
        }
        
        x <- new("Anno", annoData=aD, termData=tD, domainData=dD)
    }
)


################################################################################
################################################################################
#' @title Definition for S4 class Eoutput
#' @description \code{Eoutput} is an S4 class to store output from enrichment analysis by \code{\link{dcEnrichment}}.
#' @return Class Eoutput
#' @slot domain A character specifying the domain identity
#' @slot ontology A character specifying the ontology identity
#' @slot term_info A data.frame of nTerm X 5 containing term information, where nTerm is the number of terms in consideration, and the 5 columns are "term_id" (i.e. "Term ID"), "term_name" (i.e. "Term Name"), "namespace" (i.e. "Term Namespace"), "distance" (i.e. "Term Distance") and "IC" (i.e. "Information Content for the term based on annotation frequency by it")
#' @slot anno A list of terms, each storing annotated domains. Always, terms are identified by "term_id" and domain members identified by their ids (e.g. sunids for SCOP domains)
#' @slot data A vector containing input data in \code{\link{dcEnrichment}}. It is not always the same as the input data as only those mappable are retained
#' @slot overlap A list of terms, each storing domains overlapped between domains annotated by a term and domains in the input data (i.e. the domains of interest). Always, terms are identified by "term_id" and domain members identified by their ids (e.g. sunids for SCOP domains)
#' @slot zscore A vector of terms, containing z-scores
#' @slot pvalue A vector of terms, containing p-values
#' @slot adjp A vector of terms, containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons
#' @section Creation:
#' An object of this class can be created via: \code{new("Eoutput", domain, ontology, term_info, anno, data, overlap, zscore, pvalue, adjp)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{zscore()}: }{retrieve the slot 'zscore' in the object}
#' \item{\code{pvalue()}: }{retrieve the slot 'pvalue' in the object}
#' \item{\code{adjp()}: }{retrieve the slot 'adjp' in the object}
#' \item{\code{view()}: }{retrieve an integrated data.frame used for viewing the object}
#' \item{\code{write()}: }{write the object into a local file}
#' }
#' Standard generic methods:
#' \itemize{
#' \item{\code{str()}: }{compact display of the content in the object}
#' \item{\code{show()}: }{abbreviated display of the object}
#' }
#' @section Access:
#' Ways to access information on this class:
#' \itemize{
#' \item{\code{showClass("Eoutput")}: }{show the class definition}
#' \item{\code{showMethods(classes="Eoutput")}: }{show the method definition upon this class}
#' \item{\code{getSlots("Eoutput")}: }{get the name and class of each slot in this class}
#' \item{\code{slotNames("Eoutput")}: }{get the name of each slot in this class}
#' \item{\code{selectMethod(f, signature="Eoutput")}: }{retrieve the definition code for the method 'f' defined in this class}
#' }
#' @import methods
#' @docType class
#' @keywords S4 classes
#' @name Eoutput-class
#' @rdname Eoutput-class
#' @seealso \code{\link{Eoutput-method}}
#' @examples
#' \dontrun{
#' # 1) load SCOP.sf (as 'InfoDataFrame' object)
#' SCOP.sf <- dcRDataLoader('SCOP.sf')
#' # randomly select 20 domains
#' data <- sample(rowNames(SCOP.sf), 20)
#' 
#' # 2) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eOutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOMF")
#' eOutput
#'
#' # 3) write into the file 'Eoutput.txt' in your local directory
#' write(eOutput, file='Eoutput.txt')
#'
#' # 4) view the top 5 significant terms
#' view(eOutput, top_num=5, sortBy="pvalue", details=TRUE)
#'
#' # 4) retrieve several slots directly
#' zscore(eOutput)[1:5]
#' pvalue(eOutput)[1:5]
#' adjp(eOutput)[1:5]
#' }

#' @rdname Eoutput-class
#' @aliases Eoutput
#' @exportClass Eoutput
setClass(
    Class="Eoutput",
    representation(
        domain   = "character",
        ontology = "character",
        term_info= "data.frame",
        anno     = "list",
        data     = "vector",
        overlap  = "vector",
        zscore   = "vector",
        pvalue   = "vector",
        adjp     = "vector"
    ),
    prototype = prototype(
        domain   = "domain",
        ontology = "ontology",
        term_info= data.frame(),
        anno     = list(),
        data     = vector(),
        overlap  = vector(),
        zscore   = vector(),
        pvalue   = vector(),
        adjp     = vector()
    ),
    validity = function(object){
        if(!is.data.frame(object@term_info)){
            return("term_info is not data.frame")
        }else{
            return(TRUE)
        }
    }
)

########################################
#' @title Methods defined for S4 class Eoutput
#' @description Methods defined for S4 class \code{Eoutput}.
#' @param object an object of S4 class \code{Eoutput}. Usually this is an output from \code{\link{dcEnrichment}}
#' @param x an object of S4 class \code{Eoutput}. Usually this is an output from \code{\link{dcEnrichment}}
#' @param top_num the maximum number (5, by default) of terms will be viewed. If NULL or NA, all terms will be viewed (this can be used for the subsequent saving)
#' @param sortBy which statistics will be used for sorting and viewing terms. It can be "pvalue" for p value, "adjp" for adjusted p value, "zscore" for enrichment z-score, "nAnno" for the number in domains annotated by a term, "nOverlap" for the number in overlaps, and "none" for ordering simply according to ID of terms
#' @param decreasing logical to indicate whether to sort in a decreasing order. If it is null, by default it will be true for "zscore", "nAnno" or "nOverlap"; otherwise false
#' @param details logical to indicate whether the detailed information of terms is also viewed. By default, it sets to TRUE for the inclusion
#' @param file a character specifying a file name written into. By default, it is 'Eoutput.txt'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return
#' view(x) returns a data frame with following components:
#' \itemize{
#'  \item{\code{term_id}: term ID}
#'  \item{\code{nAnno}: number in domains annotated by a term}
#'  \item{\code{nGroup}: number in domains from the input group}
#'  \item{\code{nOverlap}: number in overlaps}
#'  \item{\code{zscore}: enrichment z-score}
#'  \item{\code{pvalue}: p value}
#'  \item{\code{adjp}: adjusted p value}
#'  \item{\code{term_name}: term name; optional, it is only appended when "details" is true}
#'  \item{\code{term_namespace}: term namespace; optional, it is only appended when "details" is true}
#'  \item{\code{term_distance}: term distance; optional, it is only appended when "details" is true}
#' }
#' write(x) also returns the same data frame as view(x), in addition to a specified local file.
#' @docType methods
#' @keywords S4 methods
#' @name Eoutput-method
#' @rdname Eoutput-method
#' @seealso \code{\link{Eoutput-class}}

#' @rdname Eoutput-method
#' @aliases show,Eoutput-method
#' @export
setMethod("show", "Eoutput",
    function(object) {
        cat(sprintf("An object of S4 class '%s', containing following slots:", class(object)), "\n", sep="")
        cat(sprintf("  @domain: '%s'", object@domain), "\n", sep="")
        cat(sprintf("  @ontology: '%s'", object@ontology), "\n", sep="")
        cat(sprintf("  @term_info: a data.frame of %d terms X %d information", dim(object@term_info)[1],dim(object@term_info)[2]), "\n", sep="")
        cat(sprintf("  @anno: a list of %d terms, each storing annotated domains", length(object@anno)), "\n", sep="")
        cat(sprintf("  @data: a vector containing a group of %d input domains", length(object@data)), "\n", sep="")
        cat(sprintf("  @overlap: a list of %d terms, each containing domains overlapped with input domains", length(object@overlap)), "\n", sep="")
        cat(sprintf("  @zscore: a vector of %d terms, containing z-scores", length(object@zscore)), "\n", sep="")
        cat(sprintf("  @pvalue: a vector of %d terms, containing p-values", length(object@pvalue)), "\n", sep="")
        cat(sprintf("  @adjp: a vector of %d terms, containing adjusted p-values", length(object@adjp)), "\n", sep="")
        cat(sprintf("In summary, a total of %d terms ('%s') are analysed for a group of %d input domains ('%s')", dim(object@term_info)[1],object@ontology,length(object@data),object@domain), "\n", sep="")
    }
)

setGeneric("zscore", function(x) standardGeneric("zscore"))
#' @rdname Eoutput-method
#' @aliases zscore
#' @export
setMethod("zscore", "Eoutput", function(x) x@zscore)

setGeneric("pvalue", function(x) standardGeneric("pvalue"))
#' @rdname Eoutput-method
#' @aliases pvalue
#' @export
setMethod("pvalue", "Eoutput", function(x) x@pvalue)

setGeneric("adjp", function(x) standardGeneric("adjp"))
#' @rdname Eoutput-method
#' @aliases adjp
#' @export
setMethod("adjp", "Eoutput", function(x) x@adjp)

setGeneric("view", function(x, ...) standardGeneric("view"))
#' @rdname Eoutput-method
#' @aliases view
#' @export
setMethod("view", "Eoutput", 
    function(x, top_num=5, sortBy=c("pvalue","adjp","zscore","nAnno","nOverlap","none"), decreasing=NULL, details=T){
        sortBy <- match.arg(sortBy)
    
        if( is.null(top_num) || is.na(top_num) ){
            top_num <- length(x@term_info$term_id)
        }
        if ( top_num > length(x@term_info$term_id) ){
            top_num <- length(x@term_info$term_id)
        }
    
        tab <- data.frame( term_id          = x@term_info$term_id,
                           nAnno            = sapply(x@anno,length),
                           nGroup           = length(x@data),
                           nOverlap         = sapply(x@overlap,length),
                           zscore           = x@zscore,
                           pvalue           = x@pvalue,
                           adjp             = x@adjp,
                           term_name        = x@term_info$term_name,
                           term_namespace   = x@term_info$term_namespace,
                           term_distance    = x@term_info$term_distance
                          )
    
        if(details == T){
            res <- tab[,c(1:10)]
        }else{
            res <- tab[,c(1:7)]
        }
    
        if(is.null(decreasing)){
            if(sortBy=="zscore" | sortBy=="nAnno" | sortBy=="nOverlap"){
                decreasing <- T
            }else{
                decreasing <- F
            }
        }
    
        switch(sortBy, 
            adjp={res <- res[order(res[,7], decreasing=decreasing)[1:top_num],]},
            pvalue={res <- res[order(res[,6], decreasing=decreasing)[1:top_num],]},
            zscore={res <- res[order(res[,5], decreasing=decreasing)[1:top_num],]},
            nAnno={res <- res[order(res[,2], decreasing=decreasing)[1:top_num],]},
            nOverlap={res <- res[order(res[,4], decreasing=decreasing)[1:top_num],]},
            none={res <- res[order(res[,1], decreasing=decreasing)[1:top_num],]}
        )
        
        res
    }
)

setGeneric("write", function(x, ...) standardGeneric("write"))
#' @rdname Eoutput-method
#' @aliases write
#' @export
setMethod("write", "Eoutput", 
    function(x, file="Eoutput.txt", verbose=T){
        if(file=='' || is.na(file) || is.null(file)){
            file <- "Eoutput.txt"
        }
        
        out <- view(x, top_num=NULL, sortBy="pvalue", details=TRUE)
        
        write.table(out, file=file, col.names=T, row.names=F, sep="\t")
        
        if(verbose){
            message(sprintf("A file ('%s') has been written into your local directory ('%s')", file, getwd()), appendLF=T)
        }
        
        invisible(out)
    }
)