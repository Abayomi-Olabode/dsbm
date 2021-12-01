library(igraph)
library(dynsbm)
#library(RJSONIO)
library(jsonlite)

#####################################################
# params to change
#setwd("/home/abayomi/Abayomi/skycross-recomb/tn93/")
#replicate <- 3
#####################################################

#start time
# start.time <- Sys.time()

# run analysis for all replicates
setwd('/home/abayomi/Abayomi/skycross-recomb/new_protocol/')

for (i in seq(1, 3, 1)) {
  
  for (j in seq(1, 100, 1)) {
    #files <- list.files(paste0("./", i, "bps/", j, "/mafft/"))
    
    files <- Sys.glob(paste0("./", i, "bps/", j, "/tn93/", 'rstep_*.tn93.csv'))
    slices <- as.integer(sapply(files, function(x) gsub(".*sw([0-9]+).*", "\\1", x)))
    files <- files[order(slices)]
    
    
    # Read TN93 files and get the maximum distance value
    # please make sure that they are comma separated files (or modify this script to apply whichever separator)
    list.m <- lapply(files, function(f) {
      x <- read.table(f, sep=',', header=T, na.strings=c("NA", "NaN"), row.names=1)
      quantile(as.matrix(x), probs = 0.30, na.rm = TRUE)
    })
    
    list.g <- vector("list", length = length(list.m))
    for (v in 1:length(list.m)) {
      f <- read.table(files[v], sep=',', header=T, na.strings="NA", row.names=1)
      co <- (as.numeric(list.m[[v]]))
      list.g[[v]] <- graph_from_adjacency_matrix(f <= co, weighted=TRUE, diag=FALSE, mode='undirected')
    }  

    # Print graphs
    #for (g in 1:length(list.g)) {
    # plot(list.g[[1]], vertex.size = 3, vertex.label = NA)
    #}
    
    cat(paste0("Graph Analsis Completed for Replicate  ", j, "\n\n"))
    
    # converting to R array
    all.vertices.name <- unique(sort(unlist(lapply(list.g, function(g) V(g)$name))))
    T <- length(list.g)
    N <- length(all.vertices.name)
    Y <- array(dim=c(T,N,N))
    
    
    for (t in 1:T){
      # random graph to be generated
      g.tmp <- graph.empty(n=N, directed=F)
      V(g.tmp)$name <- all.vertices.name
      Y[t,,] <- as.matrix(get.adjacency(union(g.tmp, list.g[[t]])))
    }
    
    
    # helper function for retrieving the highest ICL
    compute.icl <- function(dynsbm){    
      T <- ncol(dynsbm$membership)
      Q <- nrow(dynsbm$trans)
      N <- nrow(dynsbm$membership)
      pen <- 0.5*Q*log(N*(N-1)*T/2) + 0.25*Q*(Q-1)*T*log(N*(N-1)/2) # binary case
      if ("sigma" %in% names(dynsbm)) pen <- 2*pen # continuous case
      return(dynsbm$loglikelihood - ifelse(T>1,0.5*Q*(Q-1)*log(N*(T-1)),0) - pen)    
    }
    
    cat(paste0("Dynsbm Analysis Started for Replicate  ", j, "\n\n"))
    
    # print("Got here ........................................................ list.dynsbm")
    
    res <- tryCatch(
      {
        
        # dynSBM analysis
        list1.dynsbm <- select.dynsbm(Y, Qmin=1, Qmax=8, edge.type="binary", directed=FALSE, self.loop=FALSE,
                                      nb.cores=20, iter.max=20, nstart=25, perturbation.rate=0.2, fixed.param=FALSE, 
                                      bipartition=NULL, plot=TRUE)
        #save(list1.dynsbm, file = "dynsbm.Rdata")
        
        # manual ICL and loglikelihood plot
        Qmin <- ncol(list1.dynsbm[[1]]$trans)
        Qmax <- ncol(list1.dynsbm[[length(list1.dynsbm)]]$trans)
        
        # create output directory if it does not exist
        output.dir.prefix <- paste0("./", i, "bps/", j, "/dynsbm/")
        #output.dir <- paste0(j, "/")
        if (!dir.exists(file.path(output.dir.prefix))) {
          dir.create(file.path(output.dir.prefix), showWarnings = F)
        }
        
        # create json files for each group
        for (Q in Qmin:Qmax){
          json <- toJSON(list1.dynsbm[[Q]], pretty = TRUE)
          output.filepath <- paste0(output.dir.prefix, j, '_', Q, '_dynsbm.json')
          write(json, file= output.filepath)
        }
        
        # analyze the group with the highest ICL
        # dynsbm <- list1.dynsbm[[which.max(sapply(list1.dynsbm, compute.icl))]]
        dynsbm <- which.max(sapply(list1.dynsbm, compute.icl))
        
        # colnames(dynsbm$membership) <- paste0("window_", 1:ncol(dynsbm$membership))
        json.best <- toJSON(list1.dynsbm[[dynsbm]], pretty = TRUE)
        best.output.filepath <- paste0("./", i, "bps/", 'best-outputs/', paste0(i, "bps-", j), '_', dynsbm, '_dynsbm.json')
        write(json.best, file = best.output.filepath)
        
      },
      error=function(cond) {
        
        message(paste0("Warning occurred in '", j, "'.\n"))
        message("Here's the warning/error message:")
        message(paste0(cond, '\n'))
        file.duds <- paste0("./", i, "bps/", j, '-duds.txt')
        write(j, file=file.duds, append=T, sep='\n')
        return(NA)
        
      },
      warning=function(cond) {
        
        message(paste0(cond, '\n'))
        file.duds <- paste0("./", i, "bps/", j, '-duds.txt')
        write(j, file=file.duds, append=T, sep='\n')
        return(NULL)
        
      }
    )
    
    
    # stop and calculate time
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # xx <- time.taken
    # print(xx)
    
    #write(xx, file = paste0("time.txt"))
    
  }
  
  cat(paste0(i, " Breakpoint Data Completed\n\n"))
}
