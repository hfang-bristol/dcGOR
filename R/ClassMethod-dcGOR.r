require(Matrix)
################################################################################
################################################################################
#' @title Definition for class InfoDataFrame
#' @description \code{InfoDataFrame} has two slots: data and dimLabels.
#' @return Class InfoDataFrame
#' @slot data A data.frame containing terms (rows) and measured variables (columns).
#' @slot dimLabels A character descripting labels for rows and columns.
#' @section Creation:
#' An object of this class can be created via: \code{new("InfoDataFrame", data, dimLabels)}
#' @section Methods:
#' Class-specific methods:
#' \itemize{
#' \item{\code{dimLabels()}: }{retrieve labels used for display of rows and columns in the object}
#' \item{\code{dim()}: }{retrieve the dimension in the object}
#' \item{\code{nrow()}: }{retrieve number of rows in the object}
#' \item{\code{ncol()}: }{retrieve number of columns in the object}
#' \item{\code{rowNames()}: }{retrieve names of rows in the object}
#' \item{\code{colNames()}: }{retrieve names of columns in the object}
#' \item{\code{Data()}: }{retrieve the data in the object}
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
#' @keywords classes
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
#' @title Methods defined for class InfoDataFrame
#' @description Methods defined for class \code{InfoDataFrame}.
#' @param x an object of class \code{InfoDataFrame}
#' @param object an object of class \code{InfoDataFrame}
#' @param i an index
#' @param j an index
#' @param ... additional parameters
#' @docType methods
#' @keywords methods
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
        cat("An object of class ", class(object), "\n", sep="")
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
#' @title Definition for VIRTUAL class AnnoData
#' @description \code{AnnoData} is union of other classes: either matrix, or data.frame or dgCMatrix (a sparse matrix in the package Matrix). It is used as a virtual class
#' @return Class AnnoData
#' @import Matrix
#' @import methods
#' @docType class
#' @keywords classes
#' @name AnnoData-class
#' @rdname AnnoData-class
#' @seealso \code{\link{Anno-class}}

#' @rdname AnnoData-class
#' @aliases AnnoData
#' @exportClass AnnoData
setClassUnion("AnnoData", c("matrix", "data.frame", "dgCMatrix"))

################################################################################
################################################################################
#' @title Definition for class Anno
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
#' \item{\code{annoData()}: }{retrieve annoData in the object}
#' \item{\code{termData()}: }{retrieve termData (as class InfoDataFrame) in the object}
#' \item{\code{domainData()}: }{retrieve domainData (as class InfoDataFrame) in the object}
#' \item{\code{tData()}: }{retrieve termData in the object}
#' \item{\code{dData()}: }{retrieve domainData in the object}
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
#' @keywords classes
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
#' @title Methods defined for class Anno
#' @description Methods defined for class \code{Anno}.
#' @param x an object of class \code{Anno}
#' @param object an object of class \code{Anno}
#' @param i an index
#' @param j an index
#' @param ... additional parameters
#' @docType methods
#' @keywords methods
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
        cat("An object of class ", class(object), "\n", sep="")
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