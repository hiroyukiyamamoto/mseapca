pathbank2list <- function(tbl_pathbank, subject, id){

  index0 <- which(tbl_pathbank$`Pathway Subject`==subject)
  u_pathway <- unique(tbl_pathbank$`Pathway Name`[index0])
  
  pathway_all <- tbl_pathbank$`Pathway Name`
  compound_all <- tbl_pathbank[which(names(tbl_pathbank)==id)][[1]]

  M <- NULL
  for(i in 1:length(u_pathway)){
    index <- which(pathway_all==u_pathway[i])
    M[i] <- list(compound_all[index])
  }
  names(M) <- u_pathway
  
  return(M)
}

