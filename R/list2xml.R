list2xml <- function (filepath, M)
{
  N <- XML::xmlNode("file", attrs = c(path = filepath))
  for (i in 1:length(M)) {
    n <- names(M[i])
    c <- M[i][[1]]
    z = XML::xmlNode("compound_set", attrs = c(name = n))
    for (j in 1:length(c)) {
      z <- XML::append.xmlNode(z, XML::xmlNode("compound", c[j]))
    }
    N <- XML::append.xmlNode(N, z)
  }
  XML::saveXML(N,filepath)
}
