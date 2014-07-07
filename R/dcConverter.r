#' Function to convert an object between classes 'Onto' and 'igraph'
#'
#' \code{dcConverter} is supposed to convert an object between classes 'Onto' and 'igraph'.
#'
#' @param obj an object of class "Onto" or "igraph"
#' @param from a character specifying the class converted from. It can be one of "Onto" and "igraph"
#' @param to a character specifying the class converted to. It can be one of "igraph" and "Onto"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return an object of class "Onto" or "igraph"
#' @note none
#' @export
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{Onto-class}}
#' @include dcConverter.r
#' @examples
#' # 1) load onto.GOMF (as 'Onto' object)
#' on <- dcRDataLoader('onto.GOMF')
#' on
#'
#' # 2) convert the object from 'Onto' to 'igraph' class
#' ig <- dcConverter(on, from='Onto', to='igraph')
#' ig
#'
#' # 3) convert the object from 'igraph' to 'Onto' class
#' dcConverter(ig, from='igraph', to='Onto')


dcConverter <- function (obj, from=c("Onto","igraph"), to=c("igraph","Onto"), verbose=TRUE)
{
    
    from <- match.arg(from)
    to <- match.arg(to)
    
    if (class(obj) != from){
        stop(sprintf("The class of your input object '%s' is '%s', mismatched as you intended (from='%s').\n", deparse(substitute(obj)), class(obj), from))
    }
    
    if(from==to){
        warnings(sprintf("Since the class '%s' converted from is the same as the class '%s' converted to, it will return exactly what you input.\n", from, to))
        return(obj)
    }
    
    if(from=="igraph"){
        
        ## get node data frame
        data <- igraph::get.data.frame(obj, what="vertices")
        nodeI <- new("InfoDataFrame", data=data[,-1])
        ## get adjacency matrix
        if ("weight" %in% list.edge.attributes(obj)){
            adjM <- igraph::get.adjacency(obj, type="both", attr="weight", edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        }else{
            adjM <- igraph::get.adjacency(obj, type="both", attr=NULL, edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        }
        ## for Onto
        objConverted <- new("Onto", adjMatrix=adjM, nodeInfo=nodeI)
        
    }else if(from=="Onto"){
        
        ## node info
        nodes <- nInfo(obj)
        nodes <- data.frame(name=rownames(nodes), nodes)
        nodenames <- nodeNames(obj)
        
        ## adjacency matrix
        adjM <- adjMatrix(obj)
        tmp <- which(adjM!=0, arr.ind=T)
        
        ## weighted or not
        weight_flag <- T
        if(all(adjM[tmp]==1)){
            weight_flag <- F
        }
        if(weight_flag){
            relations <- data.frame(from=nodenames[tmp[,1]], to=nodenames[tmp[,2]], weight=adjM[tmp])
        }else{
            relations <- data.frame(from=nodenames[tmp[,1]], to=nodenames[tmp[,2]])
        }
        
        ## for igraph
        objConverted <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
        
    }
    
    if(verbose){
        message(sprintf("Your input object '%s' of class '%s' has been converted into an object of class '%s'.", deparse(substitute(obj)), from, to), appendLF=T)
    }
    
    return(objConverted)
}