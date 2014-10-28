#' Function to predict ontology terms given an input file containing domain architectures (including individual domains)
#'
#' \code{dcAlgoPredictMain} is supposed to predict ontology terms given an input file containing domain architectures (including individual domains).
#'
#' @param input.file an input file containing domain architectures (including individual domains). For example, a file containing UniProt ID and domain architectures for human proteins can be found in \url{http://dcgor.r-forge.r-project.org/data/Feature/hs.txt}. As seen in this example, the input file must contain the header (in the first row) and two columns: 1st column for 'SeqID' (actually these IDs can be anything), 2nd column for 'Architecture' (SCOP domain architectures, each represented as comma-separated domains). Alternatively, the input.file can be a matrix or data frame, assuming that input file has been read
#' @param output.file an output file containing predicted results. If not NULL, a tab-delimited text file will be also written out; otherwise, there is no output file (by default)
#' @param RData.HIS RData to load. This RData conveys two bits of information: 1) feature (domain) type; 2) ontology. It stores the hypergeometric scores (hscore) between features (individual domains or consecutive domain combinations) and ontology terms. The RData name tells which domain type and which ontology to use. It can be: SCOP sf domains/combinations (including "Feature2GOBP.sf", "Feature2GOMF.sf", "Feature2GOCC.sf", "Feature2HPPA.sf"), Pfam domains/combinations (including "Feature2GOBP.pfam", "Feature2GOMF.pfam", "Feature2GOCC.pfam", "Feature2HPPA.pfam"), InterPro domains (including "Feature2GOBP.interpro", "Feature2GOMF.interpro", "Feature2GOCC.interpro", "Feature2HPPA.interpro")
#' @param weight.method the method used how to weight predictions. It can be one of "none" (no weighting; by default), "copynum" for weighting copynumber of architectures, and "ic" for weighting information content (ic) of the term, "both" for weighting both copynumber and ic
#' @param merge.method the method used to merge predictions for each component feature (individual domains and their combinations derived from domain architecture). It can be one of "sum" for summing up, "max" for the maximum, and "sequential" for the sequential merging. The sequential merging is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} ranked highest hscore 
#' @param scale.method the method used to scale the predictive scores. It can be: "none" for no scaling, "linear" for being linearily scaled into the range between 0 and 1, "log" for the same as "linear" but being first log-transformed before being scaled. The scaling between 0 and 1 is done via: \eqn{\frac{S - S_{min}}{S_{max} - S_{min}}}, where \eqn{S_{min}} and \eqn{S_{max}} are the minimum and maximum values for \eqn{S}
#' @param feature.mode the mode of how to use the features thereof. It can be: "supradomains" for using all possible domain combinations (ie supradomains; including individual domains), "domains" for using individual domains only
#' @param slim.level whether only slim terms are returned. By defaut, it is NULL and all predicted terms will be reported. If it is specified as a vector containing any values from 1 to 4, then only slim terms at these levels will be reported. Here is the meaning of these values: '1' for very general terms, '2' for general terms, '3' for specific terms, and '4' for very specific terms
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.location the characters to tell the location of built-in RData files. By default, it remotely locates at "http://supfam.org/dcGOR/data" or "http://dcgor.r-forge.r-project.org/data". For the user equipped with fast internet connection, this option can be just left as default. But it is always advisable to download these files locally. Especially when the user needs to run this function many times, there is no need to ask the function to remotely download every time (also it will unnecessarily increase the runtime). For examples, these files (as a whole or part of them) can be first downloaded into your current working directory, and then set this option as: \eqn{RData.location="."}. If RData to load is already part of package itself, this parameter can be ignored (since this function will try to load it via function \code{data} first)
#' @return 
#' a term-named vector summarising the predicted scores
#' @note
#' When 'output.file' is specified, a tab-delimited text file is output, with the column names: the first two (the same as input file), 3rd for 'Term' (predicted ontology terms), 4th for 'Score' (along with predicted scores)
#' @export
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcAlgoPredict}}
#' @include dcAlgoPredictMain.r
#' @examples
#' \dontrun{
#' input.file <- "http://dcgor.r-forge.r-project.org/data/Feature/hs.txt"
#' output <- dcAlgoPredictMain(input.file, RData.HIS="Feature2GOMF.sf")
#' length(output)
#' output[1:10]
#' }

dcAlgoPredictMain <- function(input.file, output.file=NULL, RData.HIS=c("Feature2GOBP.sf","Feature2GOMF.sf","Feature2GOCC.sf","Feature2HPPA.sf","Feature2GOBP.pfam","Feature2GOMF.pfam","Feature2GOCC.pfam","Feature2HPPA.pfam","Feature2GOBP.interpro","Feature2GOMF.interpro","Feature2GOCC.interpro","Feature2HPPA.interpro"), weight.method=c("none","copynum","ic","both"), merge.method=c("sum","max","sequential"), scale.method=c("log","linear","none"), feature.mode=c("supradomains","domains"), slim.level=NULL, parallel=TRUE, multicores=NULL, verbose=T, RData.location="http://dcgor.r-forge.r-project.org/data")
{

    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    RData.HIS <- match.arg(RData.HIS)
    merge.method <- match.arg(merge.method)
    scale.method <- match.arg(scale.method)
    weight.method <- match.arg(weight.method)
    feature.mode <- match.arg(feature.mode)
    
    ## for the interpro, only 'domains' are supported
    if(length(grep("interpro", RData.HIS, perl=T))!=0){
        feature.mode <- "domains"
    }
    
    if(is.matrix(input.file) | is.data.frame(input.file)){
        if(verbose){
            message(sprintf("Load the input file ..."), appendLF=T)
        }
        input <- as.matrix(input.file)
    }else if(is.character(input.file) & input.file!=''){
        if(verbose){
            message(sprintf("Reading the input file '%s' ...", input.file), appendLF=T)
        }
        input <- as.matrix(utils::read.delim(input.file, header=T))
    }else{
        stop("The file 'input.file' must be provided!\n")
    }
    
    # determine the distinct architectures
    tmp <- base::table(input[,2])
    copynum <- as.numeric(tmp)
    names(copynum) <- names(tmp)

    if(verbose){
        message(sprintf("Predictions for %d sequences (%d distinct architectures) using '%s' RData, '%s' merge method, '%s' scale method and '%s' feature mode (%s) ...", nrow(input), length(copynum), RData.HIS, merge.method, scale.method, feature.mode, as.character(Sys.time())), appendLF=T)
    }
    pscore <- suppressMessages(dcAlgoPredict(data=names(copynum), RData.HIS=RData.HIS, merge.method=merge.method, scale.method=scale.method, feature.mode=feature.mode, slim.level=slim.level, parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location))
    
    if(verbose){
        message(sprintf("A summary in terms of ontology terms using '%s' weight method (%s)...", weight.method, as.character(Sys.time())), appendLF=T)
    }
    
    input_copynum <- as.numeric(copynum)
    output_list <- lapply(1:length(pscore), function(i){
        x <- pscore[[i]]
        if(!is.null(x)){
            tmp_df <- cbind(rep(input_copynum[i],length(x)), names(x), as.numeric(x))
            return(tmp_df)
        }
    })
    pscore_mat <- base::do.call(base::rbind, output_list)
    tscore <- base::split(x=as.numeric(pscore_mat[,3]), f=pscore_mat[,2])
    tcopynum <- base::split(x=as.numeric(pscore_mat[,1]), f=pscore_mat[,2])
    
    ## load RData containing information on 'hscore', 'ic' and 'slim'
    x <- suppressMessages(dcRDataLoader(RData=RData.HIS, verbose=verbose, RData.location=RData.location))
    ic <- x$ic
    
    # give a summary score in terms of ontology terms
    ## whether to weight according to copynumber of architectures
    if(weight.method=="copynum" | weight.method=="both"){
        sscore <- sapply(1:length(tscore), function(i) {
            base::sum(tscore[[i]]*tcopynum[[i]])
        })
    }else{
        sscore <- sapply(1:length(tscore), function(i) {
            base::sum(tscore[[i]])
        })
    }
    names(sscore) <- names(tscore)
    ## whether to weight according to ic of terms
    if(weight.method=="ic" | weight.method=="both"){
        ind <- match(names(sscore), names(ic))
        sscore <- sscore*ic[ind]
    }
    ## scale to the range from 0 to 1
    sscore <- (sscore-min(sscore))/(max(sscore)-min(sscore))
    #### force min(sscore) to be 0.0001
    sscore[sscore==0] <- 0.0001
    sscore <- signif(sscore, digits=4)
    ####
    
    if(!is.null(output.file)){
    
        if(verbose){
            message(sprintf("Preparations for output (%s)...", as.character(Sys.time())), appendLF=T)
        }
        
        # prepare the output
        tmp_seq <- input[,1]
        tmp_da <- input[,2]
        ind <- match(tmp_da, names(pscore))
        input_ps <- pscore[ind]
        output_list <- lapply(1:length(input_ps), function(i){
            x <- input_ps[[i]]
            if(!is.null(x)){
                tmp_df <- cbind(rep(tmp_seq[i],length(x)), rep(tmp_da[i],length(x)),names(x), as.numeric(x))
            }
        })
        output <- base::do.call(base::rbind, output_list)
        colnames(output) <- c(colnames(input), "Term", "Score")
        write.table(output, file=output.file, quote=F, row.names=F, sep="\t")
        
        if(file.exists(output.file)){
            message(sprintf("The predictions have been saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
        }
    }

    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    
    invisible(sscore)
}