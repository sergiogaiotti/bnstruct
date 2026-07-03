library(bnstruct)

set.seed(123)
sm_dataset <- asia()
sm_test_data_full <- raw.data(sm_dataset)
sm_test_rows <- sample(seq_len(nrow(sm_test_data_full)), min(250, nrow(sm_test_data_full)))
sm_test_data <- sm_test_data_full[sm_test_rows, , drop = FALSE]
sm_test_sizes <- node.sizes(sm_dataset)
sm_small_cols <- seq_len(min(6, ncol(sm_test_data)))
sm_small_data <- sm_test_data[, sm_small_cols, drop = FALSE]
sm_small_sizes <- sm_test_sizes[sm_small_cols]

# small synthetic 3-node chain: BDeu score equivalence makes the three optimal
# orientations tie, so the oriented DAG is order-dependent but the CPDAG is not.
sm_chain_data <- local({
    set.seed(7)
    N <- 400
    x1 <- ifelse(runif(N) < 0.7, 1L, 2L)
    x2 <- ifelse(runif(N) < 0.9, x1, 3L - x1)
    x3 <- ifelse(runif(N) < 0.9, x2, 3L - x2)
    cbind(x1, x2, x3)
})
sm_chain_sizes <- c(2L, 2L, 2L)

sm_skeleton <- function(dag) (dag | t(dag)) * 1L
sm_cpdag_key <- function(dag) {
    p <- nrow(dag)
    vs <- character(0)
    for (child in seq_len(p)) {
        pa <- which(dag[, child] == 1L)
        if (length(pa) >= 2) {
            for (a in seq_len(length(pa) - 1)) {
                for (b in seq(a + 1, length(pa))) {
                    x <- pa[a]; y <- pa[b]
                    if (dag[x, y] == 0L && dag[y, x] == 0L) {
                        vs <- c(vs, paste(min(x, y), child, max(x, y), sep = "-"))
                    }
                }
            }
        }
    }
    paste0(paste(sm_skeleton(dag), collapse = ""), "|", paste(sort(unique(vs)), collapse = ","))
}


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


test_that("oriented-DAG stability supports scalar and asymmetric bounds", {
    modes <- c("cons", "tight", "asym_cons", "asym_tight")
    results <- lapply(modes, function(mode) {
        bnstruct:::sm(sm_small_data, sm_small_sizes,
            stability = TRUE, stability.target = "dag", stability.bound = mode
        )
    })
    for (res in results) {
        expect_true(is.logical(res$stability$stable_dag))
    }
})


test_that("CPDAG stability returns a well-formed certificate", {
    res <- bnstruct:::sm(sm_chain_data, sm_chain_sizes, stability = TRUE, stability.target = "cpdag")
    st <- res$stability
    expect_equal(st$target, "cpdag")
    expect_true(is.logical(st$stable_cpdag))
    expect_true(is.logical(st$optimal_cpdag_unique))
    # the 3-node chain has one Markov equivalence class (no collider)
    expect_true(st$optimal_cpdag_unique)
    expect_true(st$n_optimal_dags >= 1)
})


test_that("CPDAG certificate is sound: certified => all LOO CPDAGs invariant", {
    res <- bnstruct:::sm(sm_chain_data, sm_chain_sizes, stability = TRUE, stability.target = "cpdag")
    if (isTRUE(res$stability$stable_cpdag)) {
        base_key <- sm_cpdag_key(bnstruct:::sm(sm_chain_data, sm_chain_sizes))
        for (i in seq_len(nrow(sm_chain_data))) {
            d <- bnstruct:::sm(sm_chain_data[-i, , drop = FALSE], sm_chain_sizes)
            expect_equal(sm_cpdag_key(d), base_key)
        }
    } else {
        succeed("certificate inconclusive on this dataset (allowed)")
    }
})


test_that("coupled parent margin certifies a strongly dependent 2-node CPDAG", {
    set.seed(11)
    N <- 300
    x1 <- ifelse(runif(N) < 0.6, 1L, 2L)
    x2 <- ifelse(runif(N) < 0.92, x1, 3L - x1)
    X <- cbind(x1, x2)
    storage.mode(X) <- "integer"
    ns <- c(2L, 2L)

    res <- bnstruct:::sm(X, ns, stability = TRUE, stability.target = "cpdag")
    st <- res$stability
    # with 2 nodes the only structural decision is the parent (edge) choice
    expect_true(st$stable_cpdag)
    expect_length(st$parent_fragile, 0)

    # soundness: certified => every leave-one-out CPDAG matches the base
    base_key <- sm_cpdag_key(bnstruct:::sm(X, ns))
    for (i in seq_len(N)) {
        expect_equal(sm_cpdag_key(bnstruct:::sm(X[-i, , drop = FALSE], ns)), base_key)
    }
})


test_that("oriented-DAG certificate is sound: certified => all LOO DAGs invariant", {
    res <- bnstruct:::sm(sm_small_data, sm_small_sizes,
        stability = TRUE, stability.target = "dag", stability.bound = "asym_tight"
    )
    if (isTRUE(res$stability$stable_dag)) {
        base_dag <- bnstruct:::sm(sm_small_data, sm_small_sizes)
        for (i in seq_len(nrow(sm_small_data))) {
            d <- bnstruct:::sm(sm_small_data[-i, , drop = FALSE], sm_small_sizes)
            expect_equal(d, base_dag)
        }
    } else {
        succeed("certificate inconclusive on this dataset (allowed)")
    }
})
