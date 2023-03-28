test_that("hamming distance works",{
  expect_equal(hammingdist("10", "10"), 0.0)
  expect_equal(hammingdist("00", "11"), 1.0)
  expect_equal(hammingdist("11", "01"), 0.5)
  expect_equal(hammingdist("11", "10"), 0.5)
  expect_equal(hammingdist("1?", "01"), 1.0)
  expect_equal(hammingdist("11", "1?"), 0.0)
  expect_equal(hammingdist("1-", "01"), 1.0)
  expect_equal(hammingdist("11", "1-"), 0.0)
})

test_that("euclidean distance works", {
  expect_equal(euclideandist(1:3, 1:3), 0.0)
  expect_equal(round(euclideandist(c(3,2), c(4,1)), 7), 1.4142136)
})
