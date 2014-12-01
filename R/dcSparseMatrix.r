#' Function to create a sparse matrix for an input file with three columns
#'
#' \code{dcSparseMatrix} is supposed to create a sparse matrix for an input file with three columns.
#'
#' @param input.file an input file containing three columns: 1st column for rows, 2nd for columns, and 3rd for numeric values. Alternatively, the input.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a list containing arguments and their default values
#' @note
#' This function is potentially useful when debugging. Because the developer does not have to specify default values for all arguments except those arguments are of interest
#' @export
#' @seealso \code{\link{dcAlgoPredictMain}}
#' @include dcSparseMatrix.r
#' @examples
#' # create a sparse matrix of 4 X 2
#' input.file <- rbind(c('R1','C1',1), c('R2','C1',1), c('R2','C2',1), c('R3','C2',1), c('R4','C1',1))
#' res <- dcSparseMatrix(input.file)
#' res
#' # get a full matrix
#' as.matrix(res)

dcSparseMatrix <- function(input.file, verbose=T)
{
    
    if(is.matrix(input.file) | is.data.frame(input.file)){
        if(verbose){
            message(sprintf("Load the input file ..."), appendLF=T)
        }
        x <- as.matrix(input.file)
    }else if(is.character(input.file) & input.file!=''){
        if(verbose){
            message(sprintf("Reading the input file '%s' ...", input.file), appendLF=T)
        }
        x <- as.matrix(utils::read.delim(input.file, header=F, sep="\t"))
    }else{
        return(NULL)
    }
    
    if(!is.null(x)){
        x_row <- sort(unique(x[,1]))
        x_col <- sort(unique(x[,2]))
        ind_row <- match(x[,1], x_row)
        ind_col <- match(x[,2], x_col)
        x.sparse <- Matrix::sparseMatrix(i=ind_row, j=ind_col, x=as.numeric(x[,3]), dims=c(length(x_row),length(x_col)))
        rownames(x.sparse) <- x_row
        colnames(x.sparse) <- x_col
        
        if(verbose){
            message(sprintf("There are %d entries, converted into a sparse matrix of %d X %d.", nrow(x), dim(x.sparse)[1], dim(x.sparse)[2]), appendLF=T)
        }
    }else{
        return(NULL)
    }
    
    invisible(x.sparse)
}
