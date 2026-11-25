test_that("formula parser keeps every supported term", {
  spec <- get_bison_model_spec((event | duration) ~ dyad(a, b) + age1 + sex + body.mass + (1 | group))
  expect_equal(spec$fixed, c("age1", "sex", "body.mass"))
  expect_equal(spec$random, "group")
  expect_false(spec$intercept)
  expect_equal(c(spec$node_1_name, spec$node_2_name), c("a", "b"))
})

test_that("terms ending in 1 are no longer mistaken for an intercept (issue #80)", {
  spec <- get_bison_model_spec((event | duration) ~ dyad(a, b) + var1 + x01)
  expect_equal(spec$fixed, c("var1", "x01"))
  expect_false(spec$intercept)
})

test_that("intercepts are parsed", {
  expect_true(get_bison_model_spec((event | duration) ~ 1 + dyad(a, b))$intercept)
  expect_false(get_bison_model_spec((event | duration) ~ 0 + dyad(a, b))$intercept)
  expect_true(get_bison_model_spec((event | duration) ~ 1)$intercept)
})

test_that("effort can be a column with digits in its name or a constant", {
  expect_equal(get_bison_model_spec((event | obs2) ~ dyad(a, b))$divisor_var_name, "obs2")
  spec <- get_bison_model_spec((event | 1) ~ dyad(a, b))
  expect_null(spec$divisor_var_name)
  expect_equal(spec$divisor_value, 1)
})

test_that("malformed formulas give clear errors or warnings", {
  expect_error(get_bison_model_spec(~ dyad(a, b)), "two-sided")
  expect_error(get_bison_model_spec(event ~ dyad(a, b)), "\\(event \\| effort\\)")
  expect_warning(get_bison_model_spec((event | duration) ~ dyad(a, b) + log(x)), "not supported")
})

test_that("undirected A-B and B-A rows refer to the same dyad", {
  df <- data.frame(n1 = c("A", "B", "A"), n2 = c("B", "A", "C"), event = c(1, 0, 1), duration = 1)
  info <- suppressMessages(get_bison_model_data((event | duration) ~ dyad(n1, n2), df, FALSE, "binary"))
  expect_equal(info$num_dyads, 2)
  expect_equal(info$row_dyad_ids, c(1, 1, 2))
  expect_equal(info$dyad_names, c("A <-> B", "A <-> C"))
})

test_that("directed A-B and B-A rows are different dyads", {
  df <- data.frame(n1 = c("A", "B", "A"), n2 = c("B", "A", "C"), event = c(1, 0, 1), duration = 1)
  info <- get_bison_model_data((event | duration) ~ dyad(n1, n2), df, TRUE, "binary")
  expect_equal(info$num_dyads, 3)
  expect_equal(info$dyad_names, c("A -> B", "B -> A", "A -> C"))
})

test_that("factor levels are used for nodes, including unobserved ones", {
  lv <- c("A", "B", "C", "D")
  df <- data.frame(n1 = factor(c("A", "B"), levels = lv), n2 = factor(c("B", "C"), levels = lv), event = 1, duration = 1)
  info <- get_bison_model_data((event | duration) ~ dyad(n1, n2), df, FALSE, "binary")
  expect_equal(names(info$node_to_idx), lv)
  # Factors with different levels are combined
  df$n2 <- factor(c("B", "C"))
  info <- get_bison_model_data((event | duration) ~ dyad(n1, n2), df, FALSE, "binary")
  expect_equal(info$num_nodes, 4)
})

test_that("data validation catches common mistakes", {
  df <- data.frame(n1 = c("A", "A"), n2 = c("B", "C"), event = c(1, 0), duration = c(1, 1))
  f <- (event | duration) ~ dyad(n1, n2)
  expect_error(get_bison_model_data((event | effort) ~ dyad(n1, n2), df, FALSE, "binary"), "not found")
  expect_error(get_bison_model_data(f, transform(df, event = c(2, 0)), FALSE, "binary"), "can't exceed")
  expect_error(get_bison_model_data(f, transform(df, event = c(0.5, 0)), FALSE, "binary"), "whole numbers")
  expect_error(get_bison_model_data(f, transform(df, event = c(NA, 0)), FALSE, "binary"), "missing values")
  expect_error(get_bison_model_data(f, transform(df, n2 = c("A", "C")), FALSE, "binary"), "Self-loops")
  expect_error(get_bison_model_data(f, transform(df, duration = c(0.5, 1)), FALSE, "binary"), "whole numbers")
  expect_silent(get_bison_model_data(f, transform(df, duration = c(0.5, 1)), FALSE, "count"))
})

test_that("design matrices are built for fixed and random effects", {
  df <- data.frame(
    n1 = c("A", "A", "B"), n2 = c("B", "C", "C"), event = c(1, 0, 1), duration = 1,
    age = c(1, 2, 3), sex = c("F", "M", "F"), group = c("g1", "g2", "g1")
  )
  info <- get_bison_model_data((event | duration) ~ dyad(n1, n2) + age + sex + (1 | group), df, FALSE, "binary")
  expect_equal(colnames(info$model_data$design_fixed), c("fixed_age", "fixed_sexF", "fixed_sexM"))
  expect_equal(info$model_data$design_fixed[, "fixed_age"], c(1, 2, 3))
  expect_equal(info$model_data$num_random, 2)
  expect_equal(info$model_data$num_random_groups, 1)
})

test_that("get_dyad_ids matches pairs in either order for undirected networks", {
  dyad_to_idx <- rbind(c(1, 2), c(1, 3), c(2, 3))
  expect_equal(get_dyad_ids(c(2, 3, 1), c(1, 2, 3), dyad_to_idx, FALSE), c(1, 3, 2))
  expect_equal(get_dyad_ids(c(2, 1), c(1, 2), dyad_to_idx, TRUE), c(NA, 1))
})

test_that("priors are validated against the model type", {
  expect_error(validate_priors(list(edge = "normal(0, 1)"), "binary_conjugate"), "beta distribution")
  expect_error(validate_priors(list(edgee = "normal(0, 1)"), "binary"), "aren't used")
  filled <- validate_priors(list(edge = "normal(-1, 1)"), "binary")
  expect_equal(filled$edge, "normal(-1, 1)")
  expect_equal(filled$fixed, get_default_priors("binary")$fixed)
  expect_error(get_default_priors("duration"), "not supported")
})

test_that("random effects without anything setting the overall level give a warning", {
  df <- data.frame(
    n1 = c("A", "A", "B", "A"), n2 = c("B", "C", "C", "B"), event = c(1, 0, 1, 1), duration = 1,
    g = c("x", "y", "x", "y"), sex = c("F", "M", "F", "M"), age = 1:4
  )
  data_for <- function(f) get_bison_model_data(f, df, FALSE, "binary")
  expect_warning(data_for((event | duration) ~ (1 | g)), "Add an intercept")
  expect_warning(data_for((event | duration) ~ age + (1 | g)), "Add an intercept")
  expect_silent(data_for((event | duration) ~ 1 + (1 | g)))
  expect_silent(data_for((event | duration) ~ dyad(n1, n2) + (1 | g)))
  expect_silent(data_for((event | duration) ~ sex + (1 | g)))
})

test_that("old prior lists with random_mean still work, with a warning", {
  priors <- c(get_default_priors("binary"), random_mean = "normal(0, 1)")
  expect_warning(validated <- validate_priors(priors, "binary"), "no longer used")
  expect_false("random_mean" %in% names(validated))
})

test_that("undirected data with rows in both orders give a message", {
  df <- data.frame(n1 = c("A", "B", "A"), n2 = c("B", "A", "C"), event = c(1, 0, 1), effort = 1)
  expect_message(get_bison_model_data((event | effort) ~ dyad(n1, n2), df, FALSE, "binary"), "1 dyad\\(s\\) have rows in both orders")
  expect_silent(get_bison_model_data((event | effort) ~ dyad(n1, n2), df, TRUE, "binary"))
  expect_silent(get_bison_model_data((event | effort) ~ dyad(n1, n2), df[c(1, 3), ], FALSE, "binary"))
})

test_that("numeric node IDs are ordered numerically and character IDs independently of locale", {
  df <- data.frame(a = c(1, 2, 9), b = c(2, 10, 10), event = 1, effort = 1)
  info <- get_bison_model_data((event | effort) ~ dyad(a, b), df, FALSE, "binary")
  expect_equal(names(info$node_to_idx), c("1", "2", "9", "10"))

  df_chr <- data.frame(a = c("b", "B", "a"), b = c("A", "a", "b"), event = 1, effort = 1)
  info_chr <- get_bison_model_data((event | effort) ~ dyad(a, b), df_chr, FALSE, "binary")
  expect_equal(names(info_chr$node_to_idx), c("A", "B", "a", "b"))
})

test_that("factor levels are kept when only one node column is a factor", {
  df <- data.frame(n1 = factor(c("A", "B"), levels = c("A", "B", "C", "D")), n2 = c("B", "C"), event = 1, effort = 1)
  info <- get_bison_model_data((event | effort) ~ dyad(n1, n2), df, FALSE, "binary")
  expect_equal(names(info$node_to_idx), c("A", "B", "C", "D"))
})

test_that("'- 1' removes the intercept and malformed dyad terms give an error", {
  spec <- get_bison_model_spec((event | effort) ~ dyad(a, b) - 1)
  expect_equal(spec$node_1_name, "a")
  expect_false(spec$intercept)
  expect_error(get_bison_model_spec((event | effort) ~ dyad(a, b) * x), "Couldn't read the formula term")
})

test_that("prior parameters that must be positive are checked", {
  expect_error(validate_priors(list(edge = "normal(0, -1)"), "binary"), "must be positive")
  expect_error(validate_priors(list(edge = "beta(0, 1)"), "binary_conjugate"), "must be positive")
  expect_error(validate_priors(list(random_std = "half-normal(0)"), "binary"), "must be positive")
  expect_silent(validate_priors(list(edge = "normal(-3, 1)"), "binary"))
})

test_that("integer and double node IDs give the same node names", {
  expect_equal(as_node_name(c(100000L, 2L)), c("100000", "2"))
  expect_equal(as_node_name(c(1e5, 2.5, 123456789012)), c("100000", "2.5", "123456789012"))
  df <- data.frame(a = c(100000L, 100000L, 300000L), b = c(200000, 300000, 100000), event = 1, effort = 2)
  info <- suppressMessages(get_bison_model_data((event | effort) ~ dyad(a, b), df, FALSE, "binary"))
  expect_equal(names(info$node_to_idx), c("100000", "200000", "300000"))
  expect_equal(format_percent(c(0.0005, 0.05, 0.055)), c("0.05", "5", "5.5"))
})

test_that("a leading -1 removes the intercept without a warning", {
  expect_silent(spec <- get_bison_model_spec((event | effort) ~ -1 + dyad(a, b)))
  expect_false(spec$intercept)
})

test_that("node IDs from factors, text and numbers give the same names", {
  expect_equal(as_node_name(factor(c(1e5, 2e5))), c("100000", "200000"))
  expect_equal(as_node_name(c("1e+05", "A1", "007")), c("100000", "A1", "007"))
  # Text IDs that only look like scientific notation are left alone
  expect_equal(as_node_name(c("2E1", "3e4", "1E2")), c("2E1", "3e4", "1E2"))
  # Each value is formatted on its own, so other values don't change its name
  expect_equal(as_node_name(c(123456789.3, 0.000001234))[1], as_node_name(123456789.3))
  df <- data.frame(n1 = factor(c(1e5, 1e5, 3e5)), n2 = c(2e5, 3e5, 2e5), event = 1, effort = 2)
  info <- suppressMessages(get_bison_model_data((event | effort) ~ dyad(n1, n2), df, FALSE, "binary"))
  expect_equal(names(info$node_to_idx), c("100000", "200000", "300000"))
})
