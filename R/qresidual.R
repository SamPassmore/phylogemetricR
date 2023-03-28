#' Quartet Sets: Internal function
#'
#' @param taxa A vector of the taxa used
#' @param k indicates the set size. Defaults to 4. You are unlikely to need to change this
#'
#' @return The number of k size sets
#' @export
#'
#' @examples
quartet_sets = function(taxa, k = 4) {
  combinat::combn(taxa, k)
}

#' Quartet Scores
#'
#' @param data the data matrix that will be used to calculate distances
#' @param quartet the quartets that will be analysed. Can be derived from `quartet_sets`
#' @param method
#'
#' @return The quartet analysed and the delta score for that quartet
#' @export
#'
#' @examples
quartet_score = function(data, quartet, method = "hamming") {
  dij = dist.phylogemetric(data[quartet[1], ], data[quartet[2], ], method = method)
  dkl = dist.phylogemetric(data[quartet[3], ], data[quartet[4], ], method = method)
  dik = dist.phylogemetric(data[quartet[1], ], data[quartet[3], ], method = method)
  djl = dist.phylogemetric(data[quartet[2], ], data[quartet[4], ], method = method)
  dil = dist.phylogemetric(data[quartet[1], ], data[quartet[4], ], method = method)
  djk = dist.phylogemetric(data[quartet[2], ], data[quartet[3], ], method = method)

  values = sort(c(dij + dkl, dik + djl, dil + djk), decreasing = TRUE)
  m1 = values[1]
  m2 = values[2]

  qresidual = (m1 - m2)^2

  return(list(quartet = quartet,
              qresidual = qresidual))
}


summarise_taxon_scores = function(sets, quartet_distances) {
  taxa = unique(c(sets))
  n_taxa = length(taxa)
  taxa_sums = rep(0, n_taxa)
  names(taxa_sums) = taxa
  scale = get_average_distances(quartet_distances)
  for (i in seq_along(quartet_distances)) {
    numerator = taxa_sums[sets[, i]] / scale
    taxa_sums[sets[, i]] = numerator / taxa_sums[sets[, i]]
  }
  # sum of distances / number of comparisons
  taxa_sums / choose(n_taxa - 1, 3)
}


get_average_distances = function(quartet_distances){
  q_scores = sapply(quartet_distances, "[[", "qresidual")
  sum(q_scores) / length(q_scores)
}

qresiduals = function(data, taxa, method = "hamming"){

  sets = quartet_sets(taxa)

  quartet_distances = apply(sets, 2, function(s) {
    quartet_score(data, s, method = method)
  })

  qresidual_scores = sapply(quartet_distances, "[[", "qresidual")

  scale = (sum(qresidual_scores) / length(qresidual_scores))

  qresidual_scores_std = qresidual_scores / scale

  qresidual_scores_std
}

round(qresiduals(t(MATRIX), colnames(MATRIX)), 2)
