library(bnstruct)

skip_on_cran()

# Load the dataset
datafile <- system.file("extdata", "data_auto_mpg.csv", package = "bnstruct")
dataset <- read.csv(datafile)

# Define indices for discrete and continuous variables
discrete_index <- c(2, 7, 8)
continuous_index <- c(1, 3, 4, 5, 6)
n_nodes <- ncol(dataset)
n_cases <- nrow(dataset)

# Prepare node sizes
ns <- rep(0, n_nodes)
ns[discrete_index] <- apply(dataset[, discrete_index], 2, function(x) length(unique(x)))
l_card <- max(ns)

# Discretization function
discretize_equal_width <- function(cont_vars, l_card) {
  cont_vars_copy <- cont_vars
  edges_norm <- (0:l_card) / l_card
  for (i in 1:ncol(cont_vars)) {
    edges <- edges_norm * (max(cont_vars[, i]) - min(cont_vars[, i])) + min(cont_vars[, i])
    cont_vars_copy[, i] <- cut(cont_vars[, i], edges, label = FALSE, include.lowest = TRUE)
  }
  return(cont_vars_copy)
}

# Re-categorization function
recat_discrete <- function(x) {
  u <- unique(x)
  n <- seq_along(u)
  re <- matrix(rep(-1, length(x)), length(x), 1)
  for (i in seq_along(x)) {
    re[i, 1] <- n[which(x[i] == u)]
  }
  return(re)
}

# Discretize the dataset
out.disc <- discretize_equal_width(dataset[, continuous_index], l_card)
data.disc <- dataset
data.disc[, continuous_index] <- out.disc
data.disc <- as.matrix(data.disc)
data.disc[, discrete_index] <- apply(data.disc[, discrete_index], 2, recat_discrete)
ns[continuous_index] <- apply(data.disc[, continuous_index], 2, function(x) length(unique(x)))

# Set storage modes
storage.mode(ns) <- "integer"
storage.mode(data.disc) <- "integer"

# Build adjacency matrix and parents/children lists
graph <- matrix(rep(0L, 64), nrow = 8)
graph[2, 3] <- 1L
graph[3, 5] <- 1L
graph[5, 1] <- 1L
graph[c(1, 5, 3), 7] <- 1L
graph[3, 4] <- 1L
graph[4, 6] <- 1L

out.par_ch <- adj_to_parents_children(graph)
parents_list <- out.par_ch$parents
parents_list <- lapply(parents_list, function(x) as.integer(x - 1L))
children_list <- out.par_ch$children
children_list <- lapply(children_list, function(x) as.integer(x - 1L))

# Sort continuous variables in reverse topological order
cont_data_permutation <- as.integer(c(6, 4, 1, 5, 3))
data.cont <- as.matrix(dataset[, cont_data_permutation])
storage.mode(data.cont) <- "double"
storage.mode(cont_data_permutation) <- "integer"

# Run discretization
max_cycles <- 8L
cuts <- .Call(
  "bnstruct_dvbn_discretize_all",
  data.cont, data.disc,
  as.integer(n_nodes),
  as.integer(n_cases),
  as.integer(length(continuous_index)),
  cont_data_permutation - 1L,
  as.integer(ns),
  parents_list, children_list,
  as.integer(max_cycles),
  as.integer(0)
)

# Reorder cuts
cuts_ordered <- rep(list(NA), length(continuous_index))
ordering <- order(cont_data_permutation)
for (i in 1:length(continuous_index)) {
  cuts_ordered[[i]] <- c(min(data.cont[, ordering[i]]), cuts[[ordering[i]]], max(data.cont[, ordering[i]]))
}

# Define the expected cuts
correct_cuts <- list(
  c(9.0, 15.25, 17.65, 20.9, 25.65, 28.9, 46.6),
  c(68.0, 70.5, 93.5, 109.0, 159.5, 259.0, 284.5, 455.0),
  c(46.0, 71.5, 99.0, 127.0, 230.0),
  c(1613.0, 2115.0, 2480.5, 2959.5, 3657.5, 5140.0),
  c(8.0, 12.35, 13.75, 16.05, 22.85, 24.8)
)

# Test that the cuts are correct
test_that("Discretization cuts are correct", {
  expect_equal(cuts_ordered, correct_cuts)
})
