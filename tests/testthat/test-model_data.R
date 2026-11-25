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
  info <- get_bison_model_data((event | duration) ~ dyad(n1, n2), df, FALSE, "binary")
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
