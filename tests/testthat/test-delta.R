library(ape)

set.seed(92929)

#### Binary Tests ####
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

simple_expected = list(delta_score = 0.115,
                       delta_taxon_scores = setNames(c(0.14375, 0.08125, 0.1125, 0.14375, 0.09375),
                                                    LETTERS[1:5]))

test_that("binary delta scores are calcuated correctly", {
  expect_equal(delta_score(t(MATRIX), colnames(MATRIX), parallel = FALSE),
               simple_expected)
})

#### Continuous Tests ####
## Simulated data
# strong phylogenetic signal
tree = rcoal(20)
df_strong = sapply(1:4, function(x) rTraitCont(tree))

## random data
df_r = matrix(rnorm(n = 80), nrow = 20)
rownames(df_r) = letters[1:20]

## Finches Data
contdata = read.table('extdata/DarwinsFinchesTraits.txt', sep = "\t",
                      header = 1)
taxa = contdata[,1]
contdata_df = contdata[,2:6]
rownames(contdata_df) = taxa

## Expected results
results = list(
  simulated_strong = list(
    delta_score = 0.232843424,
    delta_taxon_scores = c(t14 = 0.2359, t17 = 0.1807, t12 = 0.20939, t2 = 0.25532, t4 = 0.22805,
                           t6 = 0.22888, t10 = 0.21346, t11 = 0.20354, t13 = 0.25517, t8 = 0.25152,
                           t7 = 0.26726, t5 = 0.2838, t18 = 0.26277, t19 = 0.21954, t16 = 0.20684,
                           t20 = 0.20438, t9 = 0.2512, t1 = 0.1993, t15 = 0.2191, t3 = 0.28076
    )
  ),
  random = list(
    delta_score = 0.41334662,
    delta_taxon_scores =c(a = 0.41977, b = 0.40593, c = 0.43537, d = 0.42277, e = 0.39807,
                          f = 0.38353, g = 0.41656, h = 0.43771, i = 0.40344, j = 0.4299,
                          k = 0.42302, l = 0.39955, m = 0.39495, n = 0.40951, o = 0.43251,
                          p = 0.41122, q = 0.41905, r = 0.39518, s = 0.41742, t = 0.41147
    )
  ),
  finches = list(
    delta_score = 0.22995927,
    delta_taxon_scores = c(
      `Geospiza magnirostris` = 0.17557,
      `Geospiza conirostris` = 0.18676,
      `Geospiza difficilis` = 0.23118,
      `Geospiza scandens` = 0.31214,
      `Geospiza fortis` = 0.23514,
      `Geospiza fuliginosa` = 0.19611,
      `Cactospiza pallida` = 0.27903,
      `Certhidea olivacea` = 0.17757,
      `Camarhynchus parvulus` = 0.20788,
      `Camarhynchus pauper` = 0.26999,
      `Pinaroloxias inornata` = 0.21586,
      `Platyspiza crassirostris` = 0.24647,
      `Camarhynchus psittacula` = 0.25577
    )
  )
)

test_that("contunous =delta scores are calcuated correctly", {
  expect_equal(delta_score(df_strong, rownames(df_strong), method = "euclidean", parallel = TRUE),
               results$simulated_strong)
  expect_equal(delta_score(data = df_r, taxa = rownames(df_r), method = "euclidean", parallel = FALSE),
               results$random)
  expect_equal(delta_score(data = contdata_df, taxa = taxa, method = "euclidean", parallel = FALSE),
               results$finches)
})



