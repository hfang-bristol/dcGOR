#' Function to reconstruct ancestral discrete states using maximum parsimony algorithm
#'
#' \code{dcAncestralMP} is supposed to reconstruct ancestral discrete states using a maximum parsimony-modified Fitch algorithm. In a bottom-up manner, ancestral state for an internal node is determined if a state is shared in a majority by all its children. If two or more states in a majority are equally shared, this internal node is marked as unknown tie. For those ties, they are resolved being the same as its direct parent in a top-down manner. If the tie also occurs at the root, the state at the root is set to the last state (for example, 'present' for 2 states).
#'
#' @param x a vector of discrete states in the tips. It can be an unnamed vector; in this case, assumedly it has the same order as in the tree tips. More wisely, it is a named vector, whose names can be matched to the tip labels of the tree. The names of this input vector can be more than found in the tree labels, and they should contain all those in the tree labels
#' @param phy an object of class 'phylo'
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a list of architectures, containing three components for "transition", "states" and "rp":
#' \itemize{
#'  \item{\code{transition}: a posterior transition matrix between states}
#'  \item{\code{states}: a named vector storing states (extant and ancestral states)}
#'  \item{\code{rp}: a matrix of nodes X states, storing relative probability}
#' }
#' @note
#' This maximum parsimony algorithm for ancestral discrete state reconstruction is attributable to the basic idea as described in \url{http://sysbio.oxfordjournals.org/content/20/4/406.short}
#' @export
#' @seealso \code{\link{dcAncestralMP}}
#' @include dcAncestralMP.r
#' @examples
#' # provide the tree and states in the tips
#' tree <- "((((t10:5.03,t2:5.03):2.74,(t9:4.17,t5:4.17):3.60):2.80,(t3:4.05,t7:4.05):6.53):2.32,((t6:4.38,t1:4.38):2.18,(t8:2.17,t4:2.17):4.39):6.33);"
#' phy <- ape::read.tree(text=paste(tree, collapse=""))
#' x <- c(0, rep(1,4), rep(0,5))
#'
#' # reconstruct ancestral states
#' res <- dcAncestralMP(x, phy)
#' res
#'
#' # visualise the tree with ancestral states and their conditional probability
#' Ntip <- ape::Ntip(phy)
#' Nnode <- ape::Nnode(phy)
#' color <- c("white","gray")
#' ## main tree
#' ape::plot.phylo(phy, type="p", use.edge.length=TRUE, label.offset=1, show.tip.label=TRUE, show.node.label=FALSE)
#' ## tips
#' ape::tiplabels(pch=22, bg=color[as.numeric(x)+1], cex=2, adj=1)
#' ## internal nodes
#' ### relative probability
#' ape::nodelabels(thermo=res$relative[Ntip+1:Nnode,2:1], piecol=color[2:1], cex=0.75)
#' ### ancestral states
#' ape::nodelabels(text=res$states[Ntip+1:Nnode], node=Ntip+1:Nnode, frame="none", col="red", bg="transparent", cex=0.75)
#' #ape::nodelabels(text=phy$node.label, node=Ntip+1:Nnode, frame="none", col="red", bg="transparent", cex=0.75)

dcAncestralMP <- function(x, phy, verbose=T)
{
    
    if (class(phy) != "phylo"){
        stop("The input 'phy' must belong to the class 'phylo'!")
    }
    
    Ntip <- ape::Ntip(phy)
    Nnode <- ape::Nnode(phy)
    Ntot <- Ntip+Nnode

    # In the "postorder" order, the rows are arranged so that postorder tree traversal can be done by descending along these rows
    phy <- ape::reorder.phylo(phy, "postorder")
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    
    if (Nnode != Ntip-1){
        stop("The input 'phy' is not binary and rooted!")
    }
    
    if (!is.null(names(x))) {
    
        ind <- match(names(x), phy$tip.label)
        x <- x[ind[!is.na(ind)]]
        
        if(length(x) != Ntip){
            stop(message(sprintf("The names of input 'x' do not contain all of the tip labels of the input 'phy': %d NOT FOUND!", Ntip-length(x)), appendLF=T))
        }
    }

    if(verbose){
        message(sprintf("The input tree has '%d' tips.", Ntip), appendLF=T)
    }

    if (!is.factor(x)){
        x_tmp <- factor(x)
    }else{
        x_tmp <- x
    }
    nl <- nlevels(x_tmp)
    lvls <- levels(x_tmp)
    x_tmp <- as.integer(x_tmp)
    
    ################################################################################################
    if(verbose){
        message(sprintf("First, do maximum parsimony-modified Fitch algorithm in a bottom-up manner (%s) ...", as.character(Sys.time())), appendLF=T)
    }

    ## dimension: #nodes X #states (current node)
    Cx <- matrix(NA, nrow=Ntot, ncol=nl, dimnames=list(1:Ntot, lvls))
    
    ## for tips
    Cx[cbind(1:Ntip, x_tmp)] <- 1

    ## most recent common ancestor (MRCA) for each pair of tips and nodes
    mrca_node <- ape::mrca(phy, full=T)

    ## for internal nodes (in a postordered manner)
    for (i in seq(from=1, by=2, length.out=Nnode)) {
        
        j <- i + 1L
        cur <- e1[i]
        
        all_child <- which(mrca_node[cur,]==cur, arr.ind=T)
        all_child <- setdiff(all_child, cur) # exclude self
        
        ## only tips
        if(0){
            ind <- match(1:Ntip, all_child)
            all_child <- all_child[ind[!is.na(ind)]]
        }
        
        ### do calculation
        tmp <- Cx[all_child, ]
        ttmp <- apply(tmp, 2, function(x) sum(x,na.rm=T))
        ind <- which(ttmp==max(ttmp))
        if(length(ind)==1){
            Cx[cur, ind] <- 1
        }
        
    }
    
    anc <- apply(Cx, 1, function(x){
        tmp <- lvls[which(x==1)]
        if(length(tmp)==0){
            return(NA)
        }else{
            return(tmp)
        }
    })

    ################################################################################################
    if(verbose){
        message(sprintf("Second, resolve unknown states being the same as its direct parent in a top-down manner (%s) ...", as.character(Sys.time())), appendLF=T)
    }
    
    # break ties: the tie always follows the direct parent state (in a preorder)
    anc_final <- anc
    Cx_final <- Cx
    ties <- which(is.na(anc_final))
    if(length(ties) > 0){
        if(verbose){
            message(sprintf("\tbreak %d tie(s)", length(ties)), appendLF=T)
        }
    
        for(i in 1:length(ties)){
            child_ind <- ties[i]
            if(child_ind==Ntip+1){
                # break the tie at the root
                anc_final[child_ind] <- lvls[nl]
                Cx_final[child_ind, nl] <- 1
            }else{
                child <- names(child_ind)
                parent <- e1[match(child, e2)]
                parent_ind <- match(parent, names(anc))
                anc_final[child_ind] <- anc_final[parent_ind]
                Cx_final[child_ind, match(anc_final[parent_ind],lvls)] <- 1
            }
        }
    }
    
    ####################################################################################
    if(verbose){
        ## A summary of changes between states
        p2c <- cbind(anc_final[phy$edge[,1]], anc_final[phy$edge[,2]])
        p2c_final <- p2c
        ### for the root: being present
        if(nl==2){
            if(anc_final[e1[2*Nnode-1]] == lvls[nl]){
                p2c_final <- rbind(p2c_final, lvls)
            }
        }
        all <- paste(p2c_final[,1], "->", p2c_final[,2], sep='')
        changes <- sapply(unique(all), function(x){
            sum(x==all)
        })
        changes <- sort(changes)
        msg <- paste(names(changes),changes, sep=": ", collapse="\n")
    
        message(sprintf("In summary, the number of between-state changes:\n%s", msg), appendLF=T)
    }
    
    ### Calculate relative probability
    ## dimension: #nodes X #probability (current node)
    Lx <- matrix(0, nrow=Ntot, ncol=nl, dimnames=list(1:Ntot, lvls))
    ## for tips
    Lx[cbind(1:Ntip, x_tmp)] <- 1
    ## for internal nodes (in a postordered manner)
    for (i in seq(from=1, by=2, length.out=Nnode)) {
        j <- i + 1L
        cur <- e1[i]
        all_child <- which(mrca_node[cur,]==cur, arr.ind=T)
        all_child <- setdiff(all_child, cur) # exclude self
        ## only tips
        if(0){
            ind <- match(1:Ntip, all_child)
            all_child <- all_child[ind[!is.na(ind)]]
        }
        ### do calculation
        tmp <- Cx_final[all_child, ]
        ttmp <- apply(tmp, 2, function(x) sum(x,na.rm=T))
        Lx[cur,] <- ttmp
    }
    ## relative probability
    rp <- t(apply(Lx, 1, function(x) x/sum(x)))
    
    ### Calculate posterior transition matrix
    rate <- matrix(0, nl, nl)
    ind <- matrix(match(p2c_final, lvls),ncol=2)
    for(i in 1:nrow(ind)){
        rate[ind[i,1], ind[i,2]] <- rate[ind[i,1], ind[i,2]] + 1
    }
    colnames(rate) <- rownames(rate) <- lvls
    
    res <- list()
    res$transition <- rate
    res$states <- anc_final
    res$relative <- rp
    
    invisible(res)
}
