#' Function to assign (and evaluate) arguments with default values for an input function
#'
#' \code{dcFunArgs} is supposed to assign (and evaluate) arguments with default values for an input function.
#'
#' @param fun an input function name (character string)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a list containing arguments and their default values
#' @note
#' This function is potentially useful when debugging. Because the developer does not have to specify default values for all arguments except those arguments are of interest
#' @export
#' @seealso \code{\link{dcAlgoPredictMain}}
#' @include dcFunArgs.r
#' @examples
#' fun <- "dcAlgoPredictMain"
#' dcFunArgs(fun)

dcFunArgs <- function(fun, verbose=T)
{
    
    args_list <- base::formals(fun)
    args_names <- names(args_list)
    for(i in 1:length(args_list)){
        lft <- args_names[[i]]
        rgt <- paste(base::deparse(args_list[[i]]),collapse='')
        if(rgt!=''){
            tmp <- paste(lft, '=', rgt, sep=' ')
            eval(parse(text=tmp))
        }
    }
    
    if(verbose){
        message(sprintf("In function '%s', %d arguments have been assigned with default values.", fun, length(args_names)), appendLF=T)
    }
    
    invisible(args_list)
}
