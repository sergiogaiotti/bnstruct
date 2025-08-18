test_that("dvbn discretization one variable works", {
  n_cases <- 6
  n_nodes <- 3
  data_cont <- matrix(c(1.2, 2.5, 2.1, 4.0, 3.3, 5.1), n_cases, 1)
  data_disc <- matrix(
    c(1,2,1,2,1,2,   # node1
      1,1,2,2,1,2,   # node2
      1,2,2,1,2,1),  # node3
    n_cases, n_nodes
  )

  ns <- c(2L, 2L, 2L)

  parents <- list(integer(0), 1L, 2L)
  children <- list(2L, 3L, integer(0))

  cuts <- .Call("R_dvbn_discretize_one",
                data_cont, data_disc,
                as.integer(n_nodes), as.integer(n_cases),
                as.integer(ns),
                as.integer(0),  # cont_col (0-based)
                as.integer(0),  # target index (0-based)
                parents, children,
                as.integer(0))  # approx_parents
  expect_type(cuts, "double")
  expect_true(length(cuts) >= 0)
})

test_that("dvbn discretization all variables works", {
  n_cases <- 6
  n_nodes <- 3
  n_cont <- 1
  data_cont <- matrix(c(1.2, 2.5, 2.1, 4.0, 3.3, 5.1), n_cases, 1)
  data_disc <- matrix(
    c(1,2,1,2,1,2,
      1,1,2,2,1,2,
      1,2,2,1,2,1),
    n_cases, n_nodes
  )
  ns <- c(2L,2L,2L)
  parents <- list(integer(0), 1L, 2L)
  children <- list(2L, 3L, integer(0))

  out <- .Call("R_dvbn_discretize_all",
               data_cont, data_disc,
               as.integer(n_nodes), as.integer(n_cases),
               as.integer(n_cont), as.integer(0), # cont_index = 0
               as.integer(ns),
               parents, children,
               as.integer(10),  # max cycles
               as.integer(0))   # approx
  expect_type(out, "list")
  expect_true(length(out) == 1)
})
