# Internal functions for turning a formula and a data frame into the data
# structures used by the edge weight models.

# Remove any number of enclosing parentheses from a call, e.g. ((x)) -> x.
strip_parens <- function(expr) {
  while (is.call(expr) && identical(expr[[1]], as.name("("))) {
    expr <- expr[[2]]
  }
  expr
}

# Split the right-hand side of a formula into its additive terms.
split_plus <- function(expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("+")) && length(expr) == 3) {
    return(c(split_plus(expr[[2]]), split_plus(expr[[3]])))
  }
  # "- 1" removes the intercept, as in other R formulas.
  if (is.call(expr) && identical(expr[[1]], as.name("-")) && length(expr) == 3 &&
      identical(expr[[3]], 1)) {
    return(c(split_plus(expr[[2]]), list(0)))
  }
  if (is.call(expr) && identical(expr[[1]], as.name("-")) && length(expr) == 2 &&
      identical(expr[[2]], 1)) {
    return(list(0))
  }
  list(expr)
}

# Parse a bisonR edge model formula.
#
# The formula has the form `(event | effort) ~ terms`, where `effort` is a column
# name or a number, and terms can be:
#   * `dyad(node_1, node_2)` - edge weights for each dyad
#   * `0` or `1`             - no intercept (default) or an intercept
#   * `x`                    - a fixed effect of column `x`
#   * `(1 | g)`              - a random intercept for column `g`
#
# The formula is parsed as an R expression rather than as text, so any valid
# column name works (including names containing digits or dots).
get_bison_model_spec <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop(
      "`formula` must be a two-sided formula such as ",
      "(event | effort) ~ dyad(node_1, node_2).",
      call. = FALSE
    )
  }

  model_spec <- list(
    event_var_name = NULL,
    divisor_var_name = NULL,
    divisor_value = NULL,
    node_1_name = NULL,
    node_2_name = NULL,
    intercept = FALSE,
    fixed = character(0),
    random = character(0)
  )

  # Left hand side: (event | effort)
  lhs <- strip_parens(formula[[2]])
  if (!is.call(lhs) || !identical(lhs[[1]], as.name("|")) || length(lhs) != 3) {
    stop(
      "The left-hand side of the formula must have the form (event | effort), ",
      "for example (event | effort) ~ dyad(node_1, node_2).",
      call. = FALSE
    )
  }
  event_expr <- strip_parens(lhs[[2]])
  divisor_expr <- strip_parens(lhs[[3]])
  if (!is.name(event_expr)) {
    stop("The event part of (event | effort) must be a column name.", call. = FALSE)
  }
  model_spec$event_var_name <- as.character(event_expr)
  if (is.numeric(divisor_expr) && length(divisor_expr) == 1) {
    model_spec$divisor_value <- as.numeric(divisor_expr)
  } else if (is.name(divisor_expr)) {
    model_spec$divisor_var_name <- as.character(divisor_expr)
  } else {
    stop("The effort part of (event | effort) must be a column name or a number.", call. = FALSE)
  }

  # Right hand side
  for (term in split_plus(formula[[3]])) {
    inner <- strip_parens(term)
    if (is.numeric(term) && length(term) == 1 && term %in% c(0, 1)) {
      model_spec$intercept <- (term == 1)
    } else if (is.name(term)) {
      model_spec$fixed <- c(model_spec$fixed, as.character(term))
    } else if (is.call(term) && identical(term[[1]], as.name("dyad"))) {
      if (length(term) != 3 || !is.name(term[[2]]) || !is.name(term[[3]])) {
        stop("dyad() must be given exactly two column names, e.g. dyad(node_1, node_2).", call. = FALSE)
      }
      model_spec$node_1_name <- as.character(term[[2]])
      model_spec$node_2_name <- as.character(term[[3]])
    } else if (is.call(inner) && identical(inner[[1]], as.name("|")) && length(inner) == 3 &&
               identical(inner[[2]], 1) && is.name(inner[[3]])) {
      model_spec$random <- c(model_spec$random, as.character(inner[[3]]))
    } else if (grepl("dyad(", deparse1(term), fixed = TRUE)) {
      stop(
        "Couldn't read the formula term \"", deparse1(term), "\". Terms are combined with `+`, ",
        "e.g. (event | effort) ~ dyad(node_1, node_2) + x.",
        call. = FALSE
      )
    } else {
      warning(
        "Formula term \"", deparse1(term), "\" not supported by bisonR and will be ignored. ",
        "Supported terms are dyad(node_1, node_2), 0/1, column names and (1 | group).",
        call. = FALSE
      )
    }
  }

  model_spec
}

# Convert node IDs to character names. Numbers are written in full (never in
# scientific notation) and without trailing zeros, so that integer and double
# columns holding the same IDs (e.g. 100000L and 1e5) give the same name.
as_node_name <- function(x) {
  if (is.numeric(x)) {
    # Format each distinct value on its own: format() on a whole vector uses a
    # common number of digits, which could give one ID different names.
    values <- unique(x[!is.na(x)])
    names <- vapply(values, format, character(1), scientific = FALSE, trim = TRUE,
                    digits = 15, drop0trailing = TRUE)
    return(unname(names[match(x, values)]))
  }
  x <- as.character(x)
  # Labels that R itself wrote in scientific notation (e.g. factor(1e5) has the
  # label "1e+05") are written in full. Only R's own format is converted, so
  # IDs such as "2E1" or "3e4" are left as they are.
  scientific <- !is.na(x) & grepl("^-?[0-9](\\.[0-9]+)?e[+-][0-9]+$", x)
  if (any(scientific)) {
    x[scientific] <- as_node_name(as.numeric(x[scientific]))
  }
  x
}

# Sort node names: numerically if every name is a number (so 2 comes before 10),
# otherwise alphabetically in a locale-independent (C) order, so results don't
# depend on the computer's language settings.
sort_node_names <- function(names) {
  as_numbers <- suppressWarnings(as.numeric(names))
  if (!anyNA(as_numbers)) {
    return(names[order(as_numbers)])
  }
  sort(names, method = "radix")
}

# Map pairs of node indices to dyad IDs.
#
# `dyad_to_idx` is a two-column matrix of node indices, one row per dyad. For
# undirected networks the order of the two nodes doesn't matter. Returns NA for
# pairs that aren't in `dyad_to_idx`.
get_dyad_ids <- function(node_id_1, node_id_2, dyad_to_idx, directed) {
  node_id_1 <- as.integer(node_id_1)
  node_id_2 <- as.integer(node_id_2)
  dyad_to_idx <- matrix(as.integer(dyad_to_idx), ncol = 2)
  if (directed) {
    query <- paste(node_id_1, node_id_2, sep = "_")
    reference <- paste(dyad_to_idx[, 1], dyad_to_idx[, 2], sep = "_")
  } else {
    query <- paste(pmin(node_id_1, node_id_2), pmax(node_id_1, node_id_2), sep = "_")
    reference <- paste(
      pmin(dyad_to_idx[, 1], dyad_to_idx[, 2]),
      pmax(dyad_to_idx[, 1], dyad_to_idx[, 2]),
      sep = "_"
    )
  }
  match(query, reference)
}

# Check that the columns used in a model exist and contain sensible values.
validate_bison_data <- function(model_spec, observations, model_type) {
  if (!is.data.frame(observations)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (nrow(observations) == 0) {
    stop("`data` has no rows.", call. = FALSE)
  }

  needed <- c(
    model_spec$event_var_name, model_spec$divisor_var_name,
    model_spec$node_1_name, model_spec$node_2_name,
    model_spec$fixed, model_spec$random
  )
  missing_cols <- setdiff(needed, colnames(observations))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  for (col in needed) {
    if (anyNA(observations[[col]])) {
      stop("Column `", col, "` contains missing values (NA). Remove or impute them first.", call. = FALSE)
    }
  }

  event <- observations[[model_spec$event_var_name]]
  if (!is.numeric(event) || any(event < 0) || any(event != round(event))) {
    stop(
      "Column `", model_spec$event_var_name, "` must contain non-negative whole numbers ",
      "(0/1 or counts of events). For durations, see convert_duration_to_binary().",
      call. = FALSE
    )
  }

  if (!is.null(model_spec$divisor_var_name)) {
    divisor <- observations[[model_spec$divisor_var_name]]
  } else {
    divisor <- rep(model_spec$divisor_value, nrow(observations))
  }
  if (!is.numeric(divisor) || any(divisor < 0)) {
    stop("The effort (`", model_spec$divisor_var_name, "`) must be non-negative numbers.", call. = FALSE)
  }

  if (model_type == "binary") {
    if (any(divisor != round(divisor))) {
      stop(
        "For binary models the effort (number of sampling periods) must be whole numbers. ",
        "Use a count model for continuous observation effort.",
        call. = FALSE
      )
    }
    if (any(event > divisor)) {
      bad <- which(event > divisor)[1]
      stop(
        "For binary models the number of events can't exceed the number of sampling periods, ",
        "but row ", bad, " has ", event[bad], " events and effort ", divisor[bad], ".",
        call. = FALSE
      )
    }
  }
  if (model_type == "count" && any(event > 0 & divisor == 0)) {
    stop("Some rows have events but zero observation effort.", call. = FALSE)
  }

  invisible(TRUE)
}

# Build the data used by the edge models from a formula and a data frame.
get_bison_model_data <- function(formula, observations, directed, model_type) {
  model_spec <- get_bison_model_spec(formula)
  observations <- as.data.frame(observations)
  validate_bison_data(model_spec, observations, model_type)
  num_rows <- nrow(observations)

  if (!is.null(model_spec$node_1_name)) {
    node_1_names <- observations[[model_spec$node_1_name]]
    node_2_names <- observations[[model_spec$node_2_name]]

    # Work out the set of nodes. If both columns are factors with the same
    # levels, the levels (and their order) are used, so unobserved nodes are kept.
    # Otherwise the nodes are the factor levels of either column plus the
    # observed names, sorted with sort_node_names().
    if (is.factor(node_1_names) && is.factor(node_2_names) &&
        identical(levels(node_1_names), levels(node_2_names))) {
      node_levels <- unique(as_node_name(levels(node_1_names)))
    } else {
      factor_levels <- as_node_name(c(
        if (is.factor(node_1_names)) levels(node_1_names),
        if (is.factor(node_2_names)) levels(node_2_names)
      ))
      node_levels <- sort_node_names(unique(c(
        factor_levels, as_node_name(node_1_names), as_node_name(node_2_names)
      )))
    }
    node_1_factor <- factor(as_node_name(node_1_names), levels = node_levels)
    node_2_factor <- factor(as_node_name(node_2_names), levels = node_levels)

    node_to_idx <- seq_along(node_levels)
    names(node_to_idx) <- node_levels
    num_nodes <- length(node_to_idx)

    node_1_idx <- as.integer(node_1_factor)
    node_2_idx <- as.integer(node_2_factor)

    if (any(node_1_idx == node_2_idx)) {
      bad <- which(node_1_idx == node_2_idx)[1]
      stop(
        "Row ", bad, " has the same node (", node_levels[node_1_idx[bad]], ") on both sides of the dyad. ",
        "Self-loops aren't supported; remove these rows.",
        call. = FALSE
      )
    }

    # Undirected dyads are stored with the lower node index first, so that
    # A-B and B-A rows refer to the same dyad.
    if (directed) {
      dyad_a <- node_1_idx
      dyad_b <- node_2_idx
      sep <- " -> "
    } else {
      dyad_a <- pmin(node_1_idx, node_2_idx)
      dyad_b <- pmax(node_1_idx, node_2_idx)
      sep <- " <-> "
      # Tell the user if the data list pairs in both orders, which can mean the
      # data are directed (e.g. actor and receiver).
      ordered_keys <- unique(paste(node_1_idx, node_2_idx, sep = "_"))
      reversed_keys <- paste(node_2_idx, node_1_idx, sep = "_")
      num_both <- length(intersect(ordered_keys, unique(reversed_keys))) / 2
      if (num_both > 0) {
        message(
          num_both, " dyad(s) have rows in both orders (e.g. A, B and B, A). Because ",
          "directed = FALSE, both orders are combined into one edge. If node_1 and node_2 ",
          "are an actor and a receiver, use directed = TRUE."
        )
      }
    }
    row_keys <- paste(dyad_a, dyad_b, sep = "_")
    unique_keys <- unique(row_keys)
    first_rows <- match(unique_keys, row_keys)
    dyad_to_idx <- cbind(dyad_a[first_rows], dyad_b[first_rows])
    rownames(dyad_to_idx) <- seq_len(nrow(dyad_to_idx))
    dyad_ids <- match(row_keys, unique_keys)

    dyad_names <- paste(node_levels[dyad_to_idx[, 1]], node_levels[dyad_to_idx[, 2]], sep = sep)
    num_dyads <- nrow(dyad_to_idx)
  } else {
    dyad_ids <- rep(1L, num_rows)
    node_to_idx <- NULL
    dyad_to_idx <- NULL
    dyad_names <- NULL
    num_nodes <- 0
    num_dyads <- 0
  }

  event <- observations[[model_spec$event_var_name]]

  # Design matrices for fixed and random effects.
  design_fixed <- matrix(0, num_rows, 0)
  design_random <- matrix(0, num_rows, 0)
  fixed_names <- character(0)
  random_names <- character(0)
  random_group_index <- character(0)
  has_factor_fixed <- FALSE

  if (model_spec$intercept) {
    design_fixed <- cbind(design_fixed, 1)
    fixed_names <- c(fixed_names, "intercept")
  }

  for (term_name in model_spec$fixed) {
    values <- observations[[term_name]]
    if (is.character(values) || is.logical(values)) {
      values <- factor(values)
    }
    if (is.factor(values)) {
      has_factor_fixed <- TRUE
      # One column per level, as in earlier versions of bisonR.
      for (term_level in levels(values)) {
        design_fixed <- cbind(design_fixed, 1 * (values == term_level))
        fixed_names <- c(fixed_names, paste0("fixed_", term_name, term_level))
      }
    } else {
      design_fixed <- cbind(design_fixed, as.numeric(values))
      fixed_names <- c(fixed_names, paste0("fixed_", term_name))
    }
  }

  for (term_name in model_spec$random) {
    values <- as.factor(observations[[term_name]])
    for (term_level in levels(values)) {
      design_random <- cbind(design_random, 1 * (values == term_level))
      random_names <- c(random_names, paste0("random_", term_name, term_level))
      random_group_index <- c(random_group_index, paste0("random_", term_name))
    }
  }
  colnames(design_fixed) <- fixed_names
  colnames(design_random) <- random_names

  if (!is.null(model_spec$divisor_var_name)) {
    divisor <- observations[[model_spec$divisor_var_name]]
  } else {
    divisor <- rep(model_spec$divisor_value, num_rows)
  }

  model_data <- list(
    num_rows = num_rows,
    event = as.integer(round(event)),
    divisor = divisor,
    dyad_ids = as.integer(dyad_ids),
    design_fixed = design_fixed,
    design_random = design_random,
    num_edges = num_dyads,
    num_fixed = ncol(design_fixed),
    num_random = ncol(design_random),
    num_random_groups = length(unique(random_group_index)),
    random_group_index = as.integer(factor(random_group_index, levels = unique(random_group_index)))
  )
  if (model_type == "binary") {
    model_data$divisor <- as.integer(round(divisor))
  }

  # Random effects are centred on 0, so without a dyad() term, an intercept or a
  # factor fixed effect nothing sets the overall rate of events.
  if (length(model_spec$random) > 0 && is.null(model_spec$node_1_name) &&
      !model_spec$intercept && !has_factor_fixed) {
    warning(
      "Random effects are centred on 0, and the model has no dyad() term, intercept or factor ",
      "fixed effect to set the overall rate of events. Add an intercept with `1`, e.g. ",
      "(event | effort) ~ 1 + (1 | group).",
      call. = FALSE
    )
  }

  list(
    model_spec = model_spec,
    model_data = model_data,
    node_to_idx = node_to_idx,
    dyad_to_idx = dyad_to_idx,
    dyad_names = dyad_names,
    num_nodes = num_nodes,
    num_dyads = num_dyads,
    row_dyad_ids = as.integer(dyad_ids)
  )
}
