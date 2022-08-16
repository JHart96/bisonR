#' Convert group-by-individual matrix to bisonR format
#'
#' @param gbi Group-by-individual (R x C) matrix corresponding to R groups of C individuals.
#' @param group_properties A vector of length R describing properties of the groups to be added to the dataframe.
#' @param individual_properties A vector of length C describing properties of the individuals to be added to the dataframe.
#' @param individual_constraints A vector of length C describing constraints of the individuals.
#'
#' @return A dataframe for use in bisonR, where each row corresponds to a possible social event between two individuals.
#' @export
convert_gbi_to_bison <- function(gbi, group_properties=NULL, individual_properties=NULL, individual_constraints=NULL) {
  # Define empty bisonR dataframe
  df = data.frame(node_1=numeric(), node_2=numeric(), event=numeric(), group_id=numeric())
  if (!is.null(group_properties)) {
    df$group_property = numeric()
  }
  if (!is.null(individual_properties)) {
    df$node_1_property = numeric()
    df$node_2_property = numeric()
  }

  # If there are no individual constraints, set the constraints vector to the same value.
  if (is.null(individual_constraints)) {
    individual_constraints <- rep(0, ncol(gbi))
  }

  # If there are no node names on the GBI, define them numerically.
  node_names <- colnames(gbi)
  if (is.null(names(gbi))) {
    node_names <- 1:ncol(gbi)
  }

  # For each group in the GBI matrix, calculate all possible dyadic events.
  for (row_idx in 1:nrow(gbi)) {
    for (i in 1:ncol(gbi)) {
      for (j in 1:ncol(gbi)) {
        if (i < j && individual_constraints[i] == individual_constraints[j]) {
          event = gbi[row_idx, i] * gbi[row_idx, j]
          new_row <- list(
            node_1=node_names[i],
            node_2=node_names[j],
            event=event,
            group_id=row_idx
          )
          if (!is.null(group_properties)) {
            new_row$group_property=group_properties[row_idx]
          }
          if (is.null(individual_constraints)) {
            new_row$node_1_property=individual_properties[i]
            new_row$node_2_property=individual_properties[j]
          }
          df[nrow(df) + 1, ] <- new_row
        }
      }
    }
  }

  df
}
