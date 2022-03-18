node_names <- c("Rey", "Leia", "Obi-Wan", "Luke", "C-3PO", "BB-8", "R2-D2", "D-O")
node_types <- c("Lifeform", "Lifeform", "Lifeform", "Lifeform", "Droid", "Droid", "Droid", "Droid")
location_names <- c("A", "B", "C", "D", "E", "F")

n <- 8

node_types_binary <- 1 * (node_types == "Lifeform")
node_types_binary %*% t(node_types_binary)

p <- matrix(runif(n^2, max=0.25), n, n) + 0.75 * node_types_binary %*% t(node_types_binary)
p <- p * upper.tri(p)


d <- matrix(sample(10:50, size=n^2, replace=TRUE), n, n)
d <- d * upper.tri(d)

beta_loc <- rnorm(6, 0, 1)

df <- data.frame(matrix(nrow=0, ncol=6))
colnames(df) <- c("node_1", "node_2", "type_1", "type_2", "event", "location")
for (i in 1:n) {
  for (j in 1:n) {
    if (i < j) {
      for (k in 1:d[i, j]) {
        location_id <- sample(1:6, size=1)
        # At least one of them was visible, did they associate?
        logit_p <- qlogis(p[i, j])
        logit_pn <- logit_p + beta_loc[location_id]
        df[nrow(df) + 1, ] <- c(node_names[i], node_names[j], node_types[i], node_types[j], rbinom(1, 1, plogis(logit_pn)), location_names[location_id])
      }
    }
  }
}

df$node_1 <- factor(df$node_1, levels=node_names)
df$node_2 <- factor(df$node_2, levels=node_names)
df$type_1 <- factor(df$type_1, levels=c("Lifeform", "Droid"))
df$type_2 <- factor(df$type_2, levels=c("Lifeform", "Droid"))
df$location <- factor(df$location, levels=location_names)
df$event <- as.integer(df$event)
binary_ex1 <- df
usethis::use_data(binary_ex1)
