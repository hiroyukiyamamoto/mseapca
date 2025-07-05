read_pathway <-
function (fullpath) {

# read XML file
a <- XML::xmlTreeParse(fullpath)
b <- XML::xmlRoot(a)

# metabolite set list
ID <- NaN
metabolite_set <- NaN
for (i in 1:length(b)){

# // ID
s <- b[[i]]
u <- XML::xmlToList(s)
v <- u[names(u)=="compound"]
id <- as.character(unlist(v))

ID[i] <- list(id)
metabolite_set[i] <- as.character(XML::xmlAttrs(b[[i]]))
}

names(ID) <- metabolite_set
return(ID)
}
