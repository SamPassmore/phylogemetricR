hammingdist <- function(a, b){
  same <- 0.0
  compared <- 0.0
  for (i in 1:nchar(a)){
    if (substr(a, i, i) %in% c("?", "-") || substr(b, i, i) %in% c("?", "-")) {
      next
    } else if (substr(a, i, i) == substr(b, i, i)) {
      same <- same + 1.0
    }
    compared <- compared + 1.0
  }
  return(1.0 - (same / compared))
}

euclideandist <- function(a, b) sqrt(sum((a - b)^2))

dist.phylogemetric = function(a, b, method = "hamming"){
  if(method == "hamming"){
    a = paste(a, collapse = "")
    b = paste(b, collapse = "")
    dd = hammingdist(a,b)
  }

  if(method == "euclidean")
    dd = euclideandist(a, b)

  return(dd)
}
