rm(list=ls(all=TRUE))

# sample data
X <- matrix(c(8.0, 3.5, 6.5, 7.5, 2.0,
              10.0, 4.0, 9.0, 6.0, 3.5),
            nrow = 5,
            dimnames = list(c("Metabolite1", "Metabolite2", "Metabolite3", "Metabolite4", "Metabolite5"),
                            c("Sample1", "Sample2")))

metabolite_set <- list(c("Metabolite1","Metabolite3","Metabolite4"),c("Metabolite2","Metabolite5"))
names(metabolite_set) <- c("Metabolite_Set1","Metabolite_Set2")

# ssGSEA
ES <- NULL
for(i in 1:ncol(X)){
  x <- X[,i]
  
  es <- NULL
  for(k in 1:length(metabolite_set)){
    
    index <- order(x,decreasing=TRUE)
    
    y <- sort(x,decreasing=TRUE)
  
    index_p <- which(names(y) %in% metabolite_set[[k]])
    index_n <- which(!names(y) %in% metabolite_set[[k]])
  
    r <- c(length(x):1)
  
    z <- NULL
    for(j in 1:length(index_p)){
      z[j] <- r[index_p[j]]/sum(r[index_p])
    }
  
    # PG
    p <- rep(0,length(y))
    p[index_p] <- z
    p <- cumsum(p)
  
    # PN
    q <- rep(0,length(y))
    q[index_n] <- 1/(length(y)-length(index_p))
    n <- cumsum(q)
  
    es[k] <- sum(p-n)
    p <- NULL; q <- NULL

  }
  ES <- cbind(ES,es/(max(es)-min(es)))
}

rownames(ES) <- names(metabolite_set)
colnames(ES) <- colnames(X)

