context("test-SpecimenIDing_Functions")


test_that("Testing Knn identification functions work", {

  TestRes <- KnnIDingSingleInd(X = c(1:5), K = 3, GroupMembership = c(rep('A', 2), rep('B', 3)), TieBreaker = 'Remove')
  TestRes2 <- KnnIDingSingleInd(X = c(1,4,1.1,5,6), K = 4, GroupMembership = c(rep('A', 2), rep('B', 3)), TieBreaker = 'Remove')

  expect_equal(length(TestRes), 2)
  expect_equal(TestRes$Weighted.Result, TestRes$Unweighted.Result)
  expect_equal(TestRes$Weighted.Result, TestRes2$Weighted.Result)
  expect_false(TestRes2$Weighted.Result==TestRes2$Unweighted.Result)

})
