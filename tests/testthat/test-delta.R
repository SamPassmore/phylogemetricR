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

test_that("delta scores are calcuated correctly", {
  expect_equal(delta_score(MATRIX, colnames(MATRIX)),
               simple_expected)
})

#### Continuous Tests ####
contdata = read.table('extdata/DarwinsFinchesTraits.txt', sep = "\t",
                      header = 1)

contdata_t = t(contdata[,2:6])
colnames(contdata_t) = contdata[,1]
delta_score(data = contdata_t, taxa = contdata[,1], method = "euclidean")
