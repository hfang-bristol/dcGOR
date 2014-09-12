#' Function to build an object of the S4 class InfoDataframe from an input file
#'
#' \code{dcBuildInfoDataFrame} is supposed to build an object of of the S4 class \code{\link{InfoDataFrame}}, given an input file. This input file can, for example, contain the domain information. 
#'
#' @param input.file an input file used to build the object. For example, a file containing SCOP domain superfamilies (sf) can be found in \url{http://supfam.org/dcGOR/data/SCOPsf/SCOP.sf.txt}. As seen in this example, the input file must contain the header (in the first row), and entries in the first column intend to be domain identities (and must be unique)
#' @param output.file an output file used to save the built object as an RData-formatted file. If NULL, this file will be saved into "InfoDataFrame.RData" in the current working local directory
#' @return 
#' Any use-specified variable that is given on the right side of the assigement sign '<-', which contains the built \code{dcBuildInfoDataFrame} object.
#' Also, an RData file specified in "output.file" is saved in the local directory.  
#' @note If there are no use-specified variable that is given on the right side of the assigement sign '<-', then no object will be loaded onto the working environment.
#' @export
#' @seealso \code{\link{InfoDataFrame}}
#' @include dcBuildInfoDataFrame.r
#' @examples
#' \dontrun{
#' # 1) build an "InfoDataFrame" object that contains information on SCOP domain superfamilies (sf)
#' SCOP.sf <- dcBuildInfoDataFrame(input.file="http://supfam.org/dcGOR/data/SCOPsf/SCOP.sf.txt", output.file="SCOP.sf.RData")
#' SCOP.sf
#'
#' # 2) build an "InfoDataFrame" object that contains information on SCOP domain families (fa)
#' SCOP.fa <- dcBuildInfoDataFrame(input.file="http://supfam.org/dcGOR/data/SCOPfa/SCOP.fa.txt", output.file="SCOP.fa.RData")
#' SCOP.fa
#' }

dcBuildInfoDataFrame <- function(input.file, output.file="InfoDataFrame.RData")
{
    if(is.null(input.file) | is.na(input.file)){
        stop("The input file must be provided!\n")
    }
    
    if(is.null(output.file)){
        warnings("Since the output file is not provided, the function will use the default output file 'InfoDataFrame.RData'!\n")   
        output.file <- "InfoDataFrame.RData"
    }
    
    ## read the input file
    domain_info <- read.delim(input.file, header=T)
    rownames(domain_info) <- domain_info[,1]
    domain_info <- domain_info[order(domain_info[,1]),] # order according to the first column
    
    x <- as(domain_info, "InfoDataFrame")

    # remove the RData extension 
    output.var <- gsub(".RData$", "", output.file, ignore.case=T, perl=T)
    output.var <- gsub(".RDat$", "", output.var, ignore.case=T, perl=T)
    output.var <- gsub(".RDa$", "", output.var, ignore.case=T, perl=T)
    
    do.call(assign, list(output.var, x))
    save(list=output.var, file=output.file)
    
    if(file.exists(output.file)){
        message(sprintf("An object of S4 class 'Onto' has been built and saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
    }
    
    invisible(x)
}