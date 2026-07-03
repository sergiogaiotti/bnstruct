library(bnstruct)

set.seed(123)
sm_dataset <- asia()
sm_test_data_full <- raw.data(sm_dataset)
sm_test_rows <- sample(seq_len(nrow(sm_test_data_full)), min(250, nrow(sm_test_data_full)))
sm_test_data <- sm_test_data_full[sm_test_rows, , drop = FALSE]
sm_test_sizes <- node.sizes(sm_dataset)
sm_small_cols <- seq_len(min(10, ncol(sm_test_data)))
sm_small_rows <- seq_len(min(5000, nrow(sm_test_data)))
sm_small_data <- sm_test_data[sm_small_rows, sm_small_cols, drop = FALSE]
sm_small_sizes <- sm_test_sizes[sm_small_cols]


test_that("sm returns a dag when stability is off", {
    dag <- bnstruct:::sm(sm_test_data, sm_test_sizes, scoring.func = 0)
    expect_true(is.matrix(dag))
    expect_equal(dim(dag), c(8, 8))
})


test_that("sm return.tables includes F", {
    res <- bnstruct:::sm(sm_test_data, sm_test_sizes, return.tables = TRUE)
    expect_true(is.list(res))
    expect_true("F" %in% names(res))
    expect_length(res$F, 2^ncol(sm_test_data))
})


test_that("sm stability errors for non-BDeu scores", {
    expect_error(
        bnstruct:::sm(sm_test_data, sm_test_sizes, scoring.func = 1, stability = TRUE),
        "BDeu"
    )
})


test_that("sm stability supports scalar and asymmetric bounds", {
    modes <- c("cons", "tight", "asym_cons", "asym_tight")
    results <- lapply(modes, function(mode) {
        bnstruct:::sm(sm_small_data, sm_small_sizes, stability = TRUE, stability.bound = mode)
    })

    for (res in results) {
        expect_true(is.logical(res$stability$stable_dag))
    }
})


test_that("sm stability certificate implies LOO invariance", {
    modes <- c("cons", "tight", "asym_cons", "asym_tight")
    results <- lapply(modes, function(mode) {
        bnstruct:::sm(sm_small_data, sm_small_sizes, stability = TRUE, stability.bound = mode)
    })

    expect_true(all(vapply(results, function(res) is.logical(res$stability$stable_dag), logical(1))))

    if (any(vapply(results, function(res) isTRUE(res$stability$stable_dag), logical(1)))) {
        base_dag <- bnstruct:::sm(sm_small_data, sm_small_sizes)
        loo_dags <- lapply(seq_len(nrow(sm_small_data)), function(i) {
            bnstruct:::sm(sm_small_data[-i, , drop = FALSE], sm_small_sizes)
        })

        for (res in results) {
            if (isTRUE(res$stability$stable_dag)) {
                for (dag_i in loo_dags) {
                    expect_equal(dag_i, base_dag)
                }
            }
        }
    }
})
