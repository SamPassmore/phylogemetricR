
#' Title
#'
#' @param taxa A vector of the taxa used
#' @param k indicates the set size. Defaults to 4. You are unlikely to need to change this
#'
#' @return The number of k size sets
#' @export
#'
#' @examples
quartet_sets = function(taxa, k = 4){
  combinat::combn(taxa, k)
}

#' Title
#'
#' @param data the data matrix that will be used to calculate distances
#' @param quartet the quartets that will be analysed. Can be derived from `quartet_sets`
#' @param method
#'
#' @return The quartet analysed and the delta score for that quartet
#' @export
#'
#' @examples
quartet_score = function(data, quartet, method = "hamming"){
  dij = dist.phylogemetric(data[,quartet[1]], data[,quartet[2]])
  dkl = dist.phylogemetric(data[,quartet[3]], data[,quartet[4]])
  dik = dist.phylogemetric(data[,quartet[1]], data[,quartet[3]])
  djl = dist.phylogemetric(data[,quartet[2]], data[,quartet[4]])
  dil = dist.phylogemetric(data[,quartet[1]], data[,quartet[4]])
  djk = dist.phylogemetric(data[,quartet[2]], data[,quartet[3]])

  values = sort(c(dij + dkl, dik + djl, dil + djk), decreasing = TRUE)
  m1 = values[1]
  m2 = values[2]
  m3 = values[3]

  # Delta score
  if((m1 - m3) == 0){
    delta = 0
  } else {
    delta = (m1 - m2)/(m1 - m3)
  }

  return(list(quartet = quartet,
              delta = delta,
              qresidual = (m1 - m2)^2)
         )
}

#' Title
#'
#' @param sets The quartet sets that will be analysed
#' @param quartet_distances The quartet delta scores that were output from `quartet_score`
#'
#' @return
#' @export
#'
#' @examples
summarise_taxon_scores = function(sets, quartet_distances){

  taxa = unique(c(sets))
  n_taxa = length(taxa)
  taxa_sums = rep(0, n_taxa)
  names(taxa_sums) = taxa
  for(i in seq_along(quartet_distances)){
    taxa_sums[sets[,i]] = taxa_sums[sets[,i]] + quartet_distances[[i]]$delta
  }

  # sum of distances / number of comparisons
  taxa_sums / choose(n_taxa - 1, 3)
}

#' Delta Scores
#'
#' @param data The data matrix used to calculate distance between taxa. Each taxa should be a column in the data.
#' @param taxa A vector of the taxa analysed. Should relate to the columns in data
#'
#' @return The delta score for the entire dataset, and the average delta score for each taxon.
#' @export
#'
#' @examples
delta_score = function(data, taxa, method = "hamming"){
  sets = quartet_sets(taxa)

  quartet_distances = apply(sets, 2, function(s){
    quartet_score(data, s, method = method)
  })

  ### Calculate quartet distances
  ## If pbapply is available, use to print a progress bar
  if(requireNamespace("pbapply", quietly = TRUE)) {
      cl <- makeCluster(2)
      pboptions(type = "txt")
      quartet_distances = pbapply::pbapply(sets, 2, function(s) {
        quartet_score(data, s, method = method)
      }, cl = cl)
      stopCluster(cl)
  } else {
    stop("Not using pbapply")
    # apply(sets, 2, function(s){
    #   quartet_distances = quartet_score(data, s, method = method)
    # })
  }

  delta = sum(unlist(lapply(quartet_distances, "[[", "delta"))) / ncol(sets)
  delta_taxon = summarise_taxon_scores(sets, quartet_distances)

  return(list(delta_score = delta, delta_taxon_scores = round(delta_taxon, 5)))
}



