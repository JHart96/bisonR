test_that("convert_gbi_to_bison works", {
  gbi <- matrix(rbinom(20 * 10, 1, 0.25), 20, 10)
  df <- convert_gbi_to_bison(gbi)

  # Correct length
  expect_true(nrow(df) == 20 * 10 * 9 * 0.5)

  # Correct number of nodes
  expect_true(length(unique(c(df$node_1, df$node_2))) == ncol(gbi))

  # Entries are correct
  results <- rep(0, nrow(df))
  for (i in 1:nrow(df)) {
    results[i] <- gbi[ df[i, ]$group_id, df[i, ]$node_1] * gbi[ df[i, ]$group_id, df[i, ]$node_2] == df[i, ]$event
  }
  expect_true(all(results == TRUE))
})
