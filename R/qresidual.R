
## q-residual function


MATRIX = data.frame(
  A = c('1', '1', '1', '1', '0', '0', '1', '1', '1', '0', '1', '1',
        '1', '1', '0', '0', '1', '1', '1', '0'),
  B = c('1', '1', '1', '1', '0', '0', '0', '1', '1', '1', '1', '1',
        '1', '1', '1', '0', '0', '1', '1', '1'),
  C = c('1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '1', '0',
        '0', '0', '0', '1', '0', '1', '1', '1'),
  D = c('1', '0', '0', '0', '0', '1', '0', '1', '1', '1', '1', '0',
        '0', '0', '0', '1', '0', '1', '1', '1'),
  E = c('1', '0', '0', '0', '0', '1', '0', '1', '0', '1', '1', '0',
        '0', '0', '0', '1', '1', '1', '1', '1')
)

# get quartet distances

quartet_sets = function(taxa, k = 4){
  combn(taxa, k)
}

quartet_score = function(data, quartet){
  dij = hammingdist(paste(data[,quartet[1]], collapse = ""),
                    paste(data[,quartet[2]], collapse = ""))

  dkl = hammingdist(paste(data[,quartet[3]], collapse = ""),
                    paste(data[,quartet[4]], collapse = ""))

  dik = hammingdist(paste(data[,quartet[1]], collapse = ""),
                    paste(data[,quartet[3]], collapse = ""))

  djl = hammingdist(paste(data[,quartet[2]], collapse = ""),
                    paste(data[,quartet[4]], collapse = ""))

  dil = hammingdist(paste(data[,quartet[1]], collapse = ""),
                    paste(data[,quartet[4]], collapse = ""))

  djk = hammingdist(paste(data[,quartet[2]], collapse = ""),
                    paste(data[,quartet[3]], collapse = ""))

  values = sort(c(dij + dkl, dik + djl, dil + djk), decreasing = TRUE)
  m1 = values[1]
  m2 = values[2]
  m3 = values[3]

  return(list(quartet = quartet,
              delta = (m1 - m2)/(m1 - m3),
              qresidual = (m1 - m2)^2)
         )
}

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

delta_score = function(data, taxa){
  sets = quartet_sets(taxa)

  quartet_distances = apply(sets, 2, function(s){
    quartet_score(data, s)
  })

  delta = sum(unlist(lapply(quartet_distances, "[[", "delta"))) / ncol(sets)
  delta_taxon = summarise_taxon_scores(sets, quartet_distances)

  return(delta, delta_taxon)
}



