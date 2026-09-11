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


test_that("sm stability runs under BIC and AIC (cpdag, both modes)", {
    for (func in c(1L, 2L)) {
        res_e <- bnstruct:::sm(sm_small_data, sm_small_sizes,
            scoring.func = func,
            stability = TRUE, stability.target = "cpdag", stability.mode = "exact"
        )
        expect_true(is.list(res_e$stability))
        expect_identical(res_e$stability$mode, "exact")
        res_b <- bnstruct:::sm(sm_small_data, sm_small_sizes,
            scoring.func = func,
            stability = TRUE, stability.target = "cpdag", stability.mode = "bound"
        )
        expect_true(is.list(res_b$stability))
        # soundness: a bound certificate implies exact stability
        if (isTRUE(res_b$stability$stable_cpdag)) {
            expect_true(res_e$stability$stable_cpdag)
        }
    }
})

test_that("score_loo_bounds deltas are nonnegative for AIC and BIC", {
    p <- ncol(sm_small_data)
    ifm <- matrix(1L, p, 2^p)
    for (i in 1:p) {
        for (s in 0:(2^p - 1)) {
            if (bitwAnd(s, 2^(i - 1)) > 0) ifm[i, s + 1] <- 0L
        }
    }
    X <- sm_small_data
    storage.mode(X) <- "integer"
    ns <- sm_small_sizes
    storage.mode(ns) <- "integer"
    for (func in c(1L, 2L)) {
        b <- .Call("bnstruct_score_loo_bounds",
            X, ns, ifm, 1.0, func,
            PACKAGE = "bnstruct"
        )
        expect_true(min(b$lower[ifm == 1]) >= 0)
    }
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


test_that("exact CPDAG mode matches the brute-force leave-one-out oracle", {
    set.seed(21)
    N <- 150
    p <- 3
    X <- matrix(0L, N, p)
    X[, 1] <- ifelse(runif(N) < 0.65, 1L, 2L)
    for (j in 2:p) X[, j] <- ifelse(runif(N) < 0.9, X[, j - 1], 3L - X[, j - 1])
    storage.mode(X) <- "integer"
    ns <- rep(2L, p)

    exact <- bnstruct:::sm(X, ns, stability = TRUE, stability.target = "cpdag", stability.mode = "exact")$stability
    expect_equal(exact$mode, "exact")

    base_key <- sm_cpdag_key(bnstruct:::sm(X, ns))
    brute_stable <- all(vapply(
        seq_len(N),
        function(i) sm_cpdag_key(bnstruct:::sm(X[-i, , drop = FALSE], ns)) == base_key,
        logical(1)
    ))
    # exact certificate is necessary AND sufficient: it must equal the oracle
    expect_equal(exact$stable_cpdag, brute_stable)
    expect_equal(exact$n_influential == 0, brute_stable)
})


test_that("exact mode flags a structurally influential observation", {
    # a near-tie dataset so that removing one specific row flips the skeleton
    set.seed(4)
    N <- 700
    p <- 3
    X <- matrix(0L, N, p)
    X[, 1] <- ifelse(runif(N) < 0.65, 1L, 2L)
    for (j in 2:p) X[, j] <- ifelse(runif(N) < 0.92, X[, j - 1], 3L - X[, j - 1])
    storage.mode(X) <- "integer"
    ns <- rep(2L, p)

    st <- bnstruct:::sm(X, ns, stability = TRUE, stability.target = "cpdag", stability.mode = "exact")$stability
    expect_false(st$stable_cpdag)
    expect_true(st$n_influential >= 1)
    # each influential record points at a real row and reports what changed
    obs <- vapply(st$influential, function(w) w$obs, integer(1))
    expect_true(all(obs >= 1 & obs <= N))
    # removing a flagged row really does change the CPDAG
    base_key <- sm_cpdag_key(bnstruct:::sm(X, ns))
    o <- st$influential[[1]]$obs
    expect_false(sm_cpdag_key(bnstruct:::sm(X[-o, , drop = FALSE], ns)) == base_key)
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


test_that("stability.mode = 'exact' is refused for the oriented-DAG target", {
    expect_error(
        bnstruct:::sm(sm_small_data, sm_small_sizes,
            stability = TRUE, stability.target = "dag", stability.mode = "exact"
        ),
        "only for stability.target = 'cpdag'",
        fixed = TRUE
    )
    # the bound mode still works for that target
    expect_true(is.list(bnstruct:::sm(sm_small_data, sm_small_sizes,
        stability = TRUE, stability.target = "dag", stability.mode = "bound"
    )$stability))
})


test_that("Lg/Ug carry no garbage for candidate sets containing the node", {
    p <- ncol(sm_small_data)
    X <- sm_small_data
    storage.mode(X) <- "integer"
    ns <- sm_small_sizes
    storage.mode(ns) <- "integer"
    ifm <- matrix(1L, p, 2^p)
    for (i in 1:p) {
        for (s in 0:(2^p - 1)) {
            if (bitwAnd(s, 2^(i - 1)) > 0) ifm[i, s + 1] <- 0L
        }
    }
    A <- .Call("bnstruct_all_fam_log_marg_lik", X, ns, ifm, 1, 0L, PACKAGE = "bnstruct")
    bps <- .Call("bnstruct_fbp", aflml = A, PACKAGE = "bnstruct")
    b <- .Call("bnstruct_score_loo_bounds", X, ns, ifm, 1, 0L, PACKAGE = "bnstruct")
    self <- which(bps == 0L) # candidate sets containing the node itself
    expect_true(length(self) > 0)
    for (fn in c("bnstruct_compute_lg_ug_cons", "bnstruct_compute_lg_ug_tight")) {
        g <- if (fn == "bnstruct_compute_lg_ug_cons") {
            .Call(fn, bps, b$lower, b$upper, ifm, PACKAGE = "bnstruct")
        } else {
            .Call(fn, A, bps, b$lower, b$upper, ifm, PACKAGE = "bnstruct")
        }
        expect_true(all(g$lower[self] == -Inf))
        expect_true(all(g$upper[self] == -Inf))
        # the entries the certificate actually reads are unaffected
        expect_true(all(is.finite(g$lower[bps != 0L])))
    }
})


test_that("stability.bound and stability.scope are rejected for the cpdag target", {
    expect_error(
        bnstruct:::sm(sm_small_data, sm_small_sizes,
            stability = TRUE, stability.target = "cpdag", stability.bound = "tight"
        ),
        "applies only to stability.target = 'dag'",
        fixed = TRUE
    )
    expect_error(
        bnstruct:::sm(sm_small_data, sm_small_sizes,
            stability = TRUE, stability.target = "cpdag", stability.scope = "path"
        ),
        "applies only to stability.target = 'dag'",
        fixed = TRUE
    )
    # not supplying them (the default) is fine, and the reported scope is truthful
    b <- bnstruct:::sm(sm_small_data, sm_small_sizes,
        stability = TRUE, stability.target = "cpdag", stability.mode = "bound"
    )$stability
    expect_identical(b$scope, "optimal_visited")
    expect_identical(b$bound_type, "asym_tight")
    e <- bnstruct:::sm(sm_small_data, sm_small_sizes,
        stability = TRUE, stability.target = "cpdag", stability.mode = "exact"
    )$stability
    expect_identical(e$scope, "all_observations")
    # and both remain accepted for the dag target
    d <- bnstruct:::sm(sm_small_data, sm_small_sizes,
        stability = TRUE, stability.target = "dag",
        stability.bound = "asym_tight", stability.scope = "path"
    )$stability
    expect_identical(d$bound_type, "asym_tight")
})

test_that("tied-optima enumeration visits each (subset, parents) state once", {
    # Full factorial design: every pair of variables is exactly independent in the
    # data, so the optimum is the empty graph and all 10! sink orderings tie. The
    # unmemoized recursion walked those orderings, hit call.cap and reported a
    # truncated enumeration, which made condition 1 refuse the certificate.
    p <- 10
    X <- as.matrix(expand.grid(rep(list(1:2), p)))
    storage.mode(X) <- "integer"
    tab <- bnstruct:::sm(X, rep(2L, p), return.tables = TRUE)
    expect_equal(sum(tab$dag), 0)
    en <- bnstruct:::sm_enumerate_optimal_dags(tab$aflml, tab$bps, tab$F, p)
    expect_false(en$truncated)
    expect_length(en$dags, 1)
    st <- bnstruct:::sm(X, rep(2L, p), stability = TRUE, stability.mode = "bound")$stability
    expect_false(st$enumeration_truncated)
    expect_true(st$optimal_cpdag_unique)
})
