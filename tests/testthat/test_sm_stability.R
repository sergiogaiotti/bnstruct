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

# ---------------------------------------------------------------------------------------
# condition 3 as a table, exact leave-one-out by row pattern, coupled sink condition
# ---------------------------------------------------------------------------------------

sm_na_data <- local({
    set.seed(11)
    X <- sm_small_data
    X[sample(length(X), 25)] <- NA_integer_
    X
})

# ---------------------------------------------------------------------------------------
# family counting on distinct row patterns (scores, removal bounds, fused pass)
# ---------------------------------------------------------------------------------------

sm_pattern_cases <- local({
    set.seed(5)
    Xd <- rbind(sm_small_data, sm_small_data[rep(1:4, 15), ], sm_na_data)
    list(list(X = sm_small_data, ns = sm_small_sizes, lay = NULL, mfl = NULL),
         list(X = Xd, ns = sm_small_sizes, lay = NULL, mfl = NULL),
         list(X = sm_na_data, ns = sm_small_sizes, lay = NULL, mfl = matrix(2L, 1, 1)),
         list(X = sm_na_data, ns = sm_small_sizes, lay = c(1L, 1L, 2L, 2L, 3L, 3L),
              mfl = { m <- matrix(0L, 3, 3); m[1, 3] <- -1L; m[2, 3] <- -1L; m[3, 3] <- 5L; m }),
         list(X = cbind(sample(1:3, 80, TRUE), sample(1:2, 80, TRUE), sample(1:4, 80, TRUE), sample(1:3, 80, TRUE)),
              ns = c(3L, 2L, 4L, 3L), lay = NULL, mfl = NULL))
})

sm_case_ifm <- function(cs) {
    p <- ncol(cs$X)
    lay <- if (is.null(cs$lay)) rep(1, p) else cs$lay
    mfl <- if (is.null(cs$mfl)) matrix(p - 1L, 1, 1) else cs$mfl
    ifm <- bnstruct:::impossible.family.mask(p, lay, mfl)
    storage.mode(ifm) <- "integer"
    ifm
}

test_that("family scores counted on row patterns equal the per-family scorer", {
    for (cs in sm_pattern_cases) {
        X <- cs$X; storage.mode(X) <- "integer"; ns <- as.integer(cs$ns); p <- ncol(X)
        ifm <- sm_case_ifm(cs)
        for (sf in 0:2) for (ess in c(1, 0.3)) {
            af <- .Call("bnstruct_all_fam_log_marg_lik", X, ns, ifm, ess, as.integer(sf), PACKAGE = "bnstruct")
            for (i in seq_len(p)) for (P in 0:(2^p - 1)) {
                if (!ifm[i, P + 1]) {
                    expect_identical(af[i, P + 1], -Inf)
                    next
                }
                pa <- which(bitwAnd(P, 2^(seq_len(p) - 1)) > 0) - 1L
                ref <- .Call("bnstruct_score_node", X, ns, i - 1L, as.integer(pa), as.numeric(sf), ess, PACKAGE = "bnstruct")
                expect_identical(af[i, P + 1], ref)
            }
        }
    }
})

test_that("removal-delta bounds on row patterns equal brute force, and the fused pass equals separate calls", {
    for (cs in sm_pattern_cases[c(2, 4, 5)]) {
        X <- cs$X; storage.mode(X) <- "integer"; ns <- as.integer(cs$ns); p <- ncol(X)
        ifm <- sm_case_ifm(cs)
        for (sf in c(0L, 2L)) {
            b <- .Call("bnstruct_score_loo_bounds", X, ns, ifm, 1.0, sf, PACKAGE = "bnstruct")
            fb <- .Call("bnstruct_all_fam_scores_bounds", X, ns, ifm, 1.0, sf, PACKAGE = "bnstruct")
            af <- .Call("bnstruct_all_fam_log_marg_lik", X, ns, ifm, 1.0, sf, PACKAGE = "bnstruct")
            expect_identical(fb$aflml, af)
            expect_identical(fb[c("lower", "upper", "lower_k", "upper_k")], b)
            rows <- unique(c(1, 2, which(!complete.cases(X))[1:3], nrow(X)))
            rows <- rows[!is.na(rows)]
            ok <- which(ifm == 1L, arr.ind = TRUE)
            ok <- ok[seq(1, nrow(ok), length.out = min(40, nrow(ok))), , drop = FALSE]
            for (o in rows) {
                a2 <- .Call("bnstruct_all_fam_log_marg_lik", X[-o, , drop = FALSE], ns, ifm, 1.0, sf, PACKAGE = "bnstruct")
                for (f in seq_len(nrow(ok))) {
                    i <- ok[f, 1]; P1 <- ok[f, 2]
                    delta <- a2[i, P1] - af[i, P1]
                    expect_true(delta >= b$lower[i, P1] - 1e-9 && delta <= b$upper[i, P1] + 1e-9)
                    k <- X[o, i]
                    if (!is.na(k)) expect_true(delta >= b$lower_k[i, P1, k] - 1e-9 && delta <= b$upper_k[i, P1, k] + 1e-9)
                }
            }
        }
    }
})

test_that("exact leave-one-out on row patterns matches re-learning with missing values and masks", {
    for (cs in sm_pattern_cases[c(3, 4, 5)]) {
        X <- cs$X; storage.mode(X) <- "integer"; ns <- as.integer(cs$ns); p <- ncol(X)
        ifm <- sm_case_ifm(cs)
        args <- list(ns)
        if (!is.null(cs$lay)) { args$layering <- cs$lay; args$max.fanin.layers <- cs$mfl }
        else if (!is.null(cs$mfl)) args$max.fanin <- cs$mfl[1, 1]
        tab <- do.call(bnstruct:::sm, c(list(X), args, return.tables = TRUE))
        PM <- .Call("bnstruct_loo_parents_score", X, ns, ifm, 1.0, tab$aflml, 0L, PACKAGE = "bnstruct")
        bit <- 2^(seq_len(p) - 1)
        rows <- unique(c(1, which(!complete.cases(X))[1], nrow(X)))
        for (o in rows[!is.na(rows)]) {
            t2 <- do.call(bnstruct:::sm, c(list(X[-o, , drop = FALSE]), args, return.tables = TRUE))
            score_loo <- sum(vapply(seq_len(p), function(v) t2$aflml[v, PM[o, v] + 1], 0))
            expect_equal(score_loo, t2$F[2^p], tolerance = 1e-9)
        }
    }
})

test_that("bnstruct_coupled_parent_margins matches the per-decision R witness", {
    for (case in list(list(sm_small_data, sm_small_sizes), list(sm_na_data, sm_small_sizes),
                      list(sm_chain_data, sm_chain_sizes))) {
        X <- case[[1]]; ns <- as.integer(case[[2]]); p <- ncol(X)
        storage.mode(X) <- "integer"
        ifm <- matrix(1L, p, 2^p)
        for (i in seq_len(p)) ifm[i, bitwAnd(0:(2^p - 1), 2^(i - 1)) > 0] <- 0L
        tab <- bnstruct:::sm(X, ns, return.tables = TRUE)
        b <- .Call("bnstruct_score_loo_bounds", X, ns, ifm, 1.0, 0L, PACKAGE = "bnstruct")
        cm <- .Call("bnstruct_coupled_parent_margins", tab$aflml, tab$bps, ifm, ns,
                    b$lower_k, b$upper_k, b$lower, b$upper, PACKAGE = "bnstruct")
        for (i in seq_len(p)) for (U in 0:(2^p - 1)) {
            if (bitwAnd(U, 2^(i - 1)) > 0) next
            w <- bnstruct:::sm_coupled_parent_witness(i, U, U, tab$aflml, tab$bps, ifm, ns,
                                                      b$lower_k, b$upper_k, b$lower, b$upper)
            if (is.null(w)) {
                expect_true(is.na(cm$margin[i, U + 1]))
            } else {
                expect_equal(cm$margin[i, U + 1], w$margin, tolerance = 1e-9)
                expect_identical(cm$margin[i, U + 1] <= 0, isTRUE(w$fragile))
            }
            dm <- bnstruct:::sm_decoupled_parent_margin(i, U, tab$aflml, tab$bps, ifm, b$lower, b$upper)
            if (is.na(dm)) expect_true(is.na(cm$decoupled[i, U + 1])) else
                expect_equal(cm$decoupled[i, U + 1], dm, tolerance = 1e-9)
        }
    }
})

test_that("exact leave-one-out gives one result per row pattern and matches re-learning", {
    X <- rbind(sm_small_data, sm_small_data[rep(1:5, 20), ])     # many duplicated rows
    storage.mode(X) <- "integer"; ns <- as.integer(sm_small_sizes); p <- ncol(X)
    ifm <- matrix(1L, p, 2^p)
    for (i in seq_len(p)) ifm[i, bitwAnd(0:(2^p - 1), 2^(i - 1)) > 0] <- 0L
    tab <- bnstruct:::sm(X, ns, return.tables = TRUE)
    PM <- .Call("bnstruct_loo_parents_score", X, ns, ifm, 1.0, tab$aflml, 0L, PACKAGE = "bnstruct")
    key <- apply(X, 1, paste, collapse = ",")
    for (k in unique(key[duplicated(key)])) {
        rows <- which(key == k)
        expect_true(all(PM[rows, , drop = FALSE] == matrix(PM[rows[1], ], length(rows), p, byrow = TRUE)))
    }
    perm <- sample(nrow(X))
    PM2 <- .Call("bnstruct_loo_parents_score", X[perm, ], ns, ifm, 1.0, tab$aflml, 0L, PACKAGE = "bnstruct")
    expect_true(all(PM2 == PM[perm, ]))
    bit <- 2^(seq_len(p) - 1)
    for (o in c(1, 2, nrow(X))) {
        d2 <- bnstruct:::sm(X[-o, ], ns)
        dk <- matrix(0L, p, p)
        for (v in seq_len(p)) if (PM[o, v] > 0) dk[bitwAnd(PM[o, v], bit) > 0, v] <- 1L
        expect_identical(sm_cpdag_key(dk), sm_cpdag_key(d2))
    }
})

test_that("bnstruct_coupling_rows equals the per-row minimum of delta differences", {
    X <- sm_na_data; storage.mode(X) <- "integer"; ns <- as.integer(sm_small_sizes); p <- ncol(X)
    ifm <- matrix(1L, p, 2^p)
    for (i in seq_len(p)) ifm[i, bitwAnd(0:(2^p - 1), 2^(i - 1)) > 0] <- 0L
    refs <- cbind(rep(0L, p), vapply(seq_len(p), function(i) as.integer(2^(i %% p)), integer(1)))
    crl <- .Call("bnstruct_coupling_rows", X, ns, ifm, 1.0, 0L, refs, PACKAGE = "bnstruct")
    cr <- array(NA_real_, c(p, 2^p, 2))
    for (t in 1:2) for (i in seq_len(p)) cr[i, , t] <- crl$values[i, , crl$ref_index[i, t] + 1]
    drow <- function(i, P) {
        pa <- which(bitwAnd(P, 2^(seq_len(p) - 1)) > 0); ok <- rowSums(is.na(X[, c(pa, i), drop = FALSE])) == 0
        q <- if (length(pa)) prod(ns[pa]) else 1; r <- ns[i]; cfg <- rep(1, nrow(X)); m <- 1
        for (u in pa) { cfg <- cfg + (X[, u] - 1) * m; m <- m * ns[u] }
        cell <- (cfg - 1) * r + X[, i]; d <- numeric(nrow(X))
        nij <- tabulate(cfg[ok], q); nijk <- tabulate(cell[ok], q * r)
        d[ok] <- log(1 / q + nij[cfg[ok]] - 1) - log(1 / (q * r) + nijk[cell[ok]] - 1)
        d
    }
    for (t in 1:2) for (i in seq_len(p)) {
        dR <- drow(i, refs[i, t])
        for (P in which(ifm[i, ] == 1L) - 1) {
            ref <- if (P == refs[i, t]) 0 else min(dR - drow(i, P))
            expect_equal(cr[i, P + 1, t], ref, tolerance = 1e-12)
        }
    }
})

test_that("coupled sink bound holds for every row, subset and competitor (brute force)", {
    X <- sm_chain_data[1:80, ]; storage.mode(X) <- "integer"
    X <- cbind(X, ifelse(runif(80) < 0.8, X[, 3], 3L - X[, 3]))
    ns <- rep(2L, 4); p <- 4; full <- 2^p - 1; bit <- 2^(0:(p - 1))
    ifm <- matrix(1L, p, 2^p)
    for (i in seq_len(p)) ifm[i, bitwAnd(0:full, bit[i]) > 0] <- 0L
    dv <- function(tab) { A <- matrix(NA_real_, 2^p, p)
        for (S in 1:full) for (j in which(bitwAnd(S, bit) > 0)) {
            U <- bitwXor(S, bit[j]); A[S + 1, j] <- tab$F[U + 1] + tab$aflml[j, tab$bps[j, U + 1]] }
        A }
    tab <- bnstruct:::sm(X, ns, return.tables = TRUE)
    pa <- vapply(seq_len(p), function(v) as.integer(sum(bit[tab$dag[, v] == 1L])), integer(1))
    crl <- .Call("bnstruct_coupling_rows", X, ns, ifm, 1.0, 0L, matrix(pa, p, 1), PACKAGE = "bnstruct")
    cc <- crl$values[, , 1]
    st <- tab$aflml; fin <- is.finite(st); st[fin] <- st[fin] - cc[fin]
    bpt <- .Call("bnstruct_fbp", aflml = st, PACKAGE = "bnstruct")
    At <- dv(list(F = .Call("bnstruct_fbs_with_scores", bps = bpt, aflml = st, PACKAGE = "bnstruct")$scores,
                  aflml = st, bps = bpt))
    for (o in seq_len(nrow(X))) {
        t2 <- bnstruct:::sm(X[-o, ], ns, return.tables = TRUE)
        Ao <- dv(t2)
        for (S in 1:full) {
            nodes <- which(bitwAnd(S, bit) > 0)
            if (any(bitwAnd(pa[nodes], bitwNot(S)) > 0)) next
            scoreG <- sum(tab$aflml[cbind(nodes, pa[nodes] + 1)])
            scoreG_o <- sum(t2$aflml[cbind(nodes, pa[nodes] + 1)])
            for (j in nodes) expect_gte(scoreG_o - Ao[S + 1, j] + 1e-7, scoreG - At[S + 1, j])
        }
    }
})

test_that("coupled sink certificate is sound and never weaker than the decoupled one", {
    for (case in list(list(sm_small_data, sm_small_sizes), list(sm_chain_data, sm_chain_sizes),
                      list(sm_test_data, sm_test_sizes))) {
        X <- case[[1]]; ns <- case[[2]]
        e <- bnstruct:::sm(X, ns, stability = TRUE, stability.mode = "exact")$stability
        d <- bnstruct:::sm(X, ns, stability = TRUE, stability.mode = "bound")$stability
        c <- bnstruct:::sm(X, ns, stability = TRUE, stability.mode = "bound", stability.sink = "coupled")$stability
        expect_identical(c$sink, "coupled")
        if (isTRUE(c$stable_cpdag)) expect_true(e$stable_cpdag)
        if (isTRUE(d$stable_cpdag)) expect_true(c$stable_cpdag)
    }
    expect_error(
        bnstruct:::sm(sm_small_data, sm_small_sizes, stability = TRUE, stability.target = "dag",
                      stability.sink = "coupled"),
        "only for stability.target = 'cpdag'", fixed = TRUE
    )
})

test_that("row-coupled parent margins are sound and never looser than the child-value ones", {
    for (case in list(list(sm_small_data, sm_small_sizes), list(sm_na_data, sm_small_sizes),
                      list(sm_chain_data, sm_chain_sizes), list(sm_test_data, sm_test_sizes))) {
        X <- case[[1]]; ns <- as.integer(case[[2]]); p <- ncol(X)
        storage.mode(X) <- "integer"
        ifm <- matrix(1L, p, 2^p)
        for (i in seq_len(p)) ifm[i, bitwAnd(0:(2^p - 1), 2^(i - 1)) > 0] <- 0L
        tab <- bnstruct:::sm(X, ns, return.tables = TRUE)
        b <- .Call("bnstruct_score_loo_bounds", X, ns, ifm, 1.0, 0L, PACKAGE = "bnstruct")
        cm <- .Call("bnstruct_coupled_parent_margins", tab$aflml, tab$bps, ifm, ns,
                    b$lower_k, b$upper_k, b$lower, b$upper, PACKAGE = "bnstruct")
        refs <- t(vapply(seq_len(p), function(i) as.integer(sort(unique(tab$bps[i, ][tab$bps[i, ] > 0] - 1L))[1:4]), integer(4)))
        refs[is.na(refs)] <- 0L
        crl <- .Call("bnstruct_coupling_rows", X, ns, ifm, 1.0, 0L, refs, PACKAGE = "bnstruct")
        pr <- .Call("bnstruct_parent_margins_rows", tab$aflml, tab$bps, ifm, crl$values, refs, crl$ref_index, PACKAGE = "bnstruct")
        filled <- which(!is.na(pr$margin))
        expect_true(length(filled) > 0)
        ok <- filled[!is.na(cm$margin[filled])]
        expect_true(all(pr$margin[ok] >= cm$margin[ok] - 1e-9))
        # soundness by brute force on a few removals
        for (o in c(1, 5, nrow(X))) {
            t2 <- bnstruct:::sm(X[-o, ], ns, return.tables = TRUE)
            for (pos in filled) {
                i <- (pos - 1) %% p + 1; U <- (pos - 1) %/% p; R <- tab$bps[i, U + 1]
                sub <- which(bitwAnd(0:(2^p - 1), bitwNot(U)) == 0 & ifm[i, ] == 1L & seq_len(2^p) != R)
                if (!length(sub)) next
                gap <- t2$aflml[i, R] - max(t2$aflml[i, sub])
                expect_gte(gap + 1e-7, pr$margin[pos])
            }
        }
    }
    for (case in list(list(sm_small_data, sm_small_sizes), list(sm_test_data, sm_test_sizes))) {
        e <- bnstruct:::sm(case[[1]], case[[2]], stability = TRUE, stability.mode = "exact")$stability
        r <- bnstruct:::sm(case[[1]], case[[2]], stability = TRUE, stability.mode = "bound",
                           stability.sink = "coupled", stability.parent = "rows")$stability
        k <- bnstruct:::sm(case[[1]], case[[2]], stability = TRUE, stability.mode = "bound",
                           stability.sink = "coupled")$stability
        if (isTRUE(r$stable_cpdag)) expect_true(e$stable_cpdag)
        if (isTRUE(k$stable_cpdag)) expect_true(r$stable_cpdag)
    }
})

