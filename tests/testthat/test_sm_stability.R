library(bnstruct)

set.seed(123)
sm_dataset <- asia()
sm_test_data_full <- raw.data(sm_dataset)
sm_test_rows <- sample(seq_len(nrow(sm_test_data_full)), min(200, nrow(sm_test_data_full)))
sm_test_data <- sm_test_data_full[sm_test_rows, , drop = FALSE]
sm_test_sizes <- node.sizes(sm_dataset)


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


test_that("sm stability supports cons and tight bounds", {
    res_cons <- bnstruct:::sm(sm_test_data, sm_test_sizes, stability = TRUE, stability.bound = "cons")
    res_tight <- bnstruct:::sm(sm_test_data, sm_test_sizes, stability = TRUE, stability.bound = "tight")
    expect_true(is.logical(res_cons$stability$stable_dag))
    expect_true(is.logical(res_tight$stability$stable_dag))
})


test_that("sm stability certificate implies LOO invariance", {
    res <- bnstruct:::sm(sm_test_data, sm_test_sizes, stability = TRUE, stability.bound ="tight")
    expect_true(is.logical(res$stability$stable_dag))

    if (isTRUE(res$stability$stable_dag)) {
        for (i in seq_len(nrow(sm_test_data))) {
            dag_i <- bnstruct:::sm(sm_test_data[-i, , drop = FALSE], sm_test_sizes)
            expect_equal(dag_i, res$dag)
        }
    }
})
