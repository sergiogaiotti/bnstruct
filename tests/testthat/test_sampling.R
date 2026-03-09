library(bnstruct)

# Learn a network from the asia dataset to use for sampling tests
dataset <- asia()
network <- learn.network(dataset, algo = "sm")

# Old (for-loop) implementation for comparison
sample.dataset.old <- function(x, n = 100, mar = 0) {
    if (mar < 0 || mar > 1) {
        mar <- 0
    }
    nn <- num.nodes(x)
    obs <- matrix(rep(0, nn * n), nrow = n, ncol = nn)
    for (i in 1:n) {
        obs[i, ] <- sample.row(x, mar)
    }
    BNDataset(obs, discreteness(x), variables(x), node.sizes(x))
}


# --- Correctness tests ---

test_that("sample.dataset returns correct dimensions", {
    sampled <- sample.dataset(network, n = 500)
    expect_equal(nrow(raw.data(sampled)), 500)
    expect_equal(ncol(raw.data(sampled)), num.nodes(network))
})

test_that("sampled values are within valid ranges", {
    sampled <- sample.dataset(network, n = 1000)
    rd <- raw.data(sampled)
    ns <- node.sizes(network)
    for (v in 1:num.nodes(network)) {
        vals <- rd[, v]
        expect_true(all(vals >= 1 & vals <= ns[v]),
            info = paste("Variable", v, "has out-of-range values")
        )
    }
})

test_that("sample.dataset marginal distributions match CPTs for root nodes", {
    set.seed(42)
    n <- 50000
    sampled <- sample.dataset(network, n = n)
    rd <- raw.data(sampled)
    dag.mat <- dag(network)
    cpts.list <- cpts(network)
    ns <- node.sizes(network)
    vars <- variables(network)

    # Find root nodes (no parents)
    roots <- which(colSums(dag.mat) == 0)

    for (node in roots) {
        observed.freq <- tabulate(rd[, node], nbins = ns[node]) / n
        expected.prob <- as.numeric(cpts.list[[node]])
        # With 50k samples, allow reasonable tolerance
        expect_equal(observed.freq, expected.prob,
            tolerance = 0.02,
            info = paste("Root node", vars[node], "marginal mismatch")
        )
    }
})

test_that("sample.dataset conditional distributions match CPTs", {
    set.seed(123)
    n <- 100000
    sampled <- sample.dataset(network, n = n)
    rd <- raw.data(sampled)
    dag.mat <- dag(network)
    cpts.list <- cpts(network)
    ns <- node.sizes(network)
    vars <- variables(network)

    # Find a node with exactly one parent for a simple conditional check
    n.parents <- colSums(dag.mat)
    nodes.one.parent <- which(n.parents == 1)

    if (length(nodes.one.parent) > 0) {
        node <- nodes.one.parent[1]
        parent <- which(dag.mat[, node] != 0)

        # Get CPT dimensions
        dv <- match(
            c(unlist(names(dimnames(cpts.list[[node]])))),
            c(vars)
        )

        for (pval in 1:ns[parent]) {
            mask <- which(rd[, parent] == pval)
            if (length(mask) > 100) {
                observed.freq <- tabulate(rd[mask, node], nbins = ns[node]) / length(mask)

                # Extract expected conditional from CPT
                # The CPT is indexed by (parent_value, node_value) or vice versa
                cpt <- cpts.list[[node]]
                node.dim.pos <- which(dv == node)
                parent.dim.pos <- which(dv == parent)

                # Build index list for array subsetting
                idx <- rep(list(TRUE), length(dv))
                idx[[parent.dim.pos]] <- pval
                expected.prob <- do.call(`[`, c(list(cpt), idx))

                expect_equal(observed.freq, as.numeric(expected.prob),
                    tolerance = 0.03,
                    info = paste("Conditional P(", vars[node], "|", vars[parent], "=", pval, ") mismatch")
                )
            }
        }
    }
})

test_that("sample.dataset handles mar parameter", {
    set.seed(99)
    n <- 5000
    mar <- 0.2
    sampled <- sample.dataset(network, n = n, mar = mar)
    rd <- raw.data(sampled)
    frac.missing <- sum(is.na(rd)) / length(rd)
    # Should be roughly 20% missing, allow generous tolerance
    expect_equal(frac.missing, mar, tolerance = 0.03)
})

test_that("vectorized sample.dataset matches old for-loop version statistically", {
    # Both should produce samples from the same distribution.
    # Compare marginal distributions from both methods.
    set.seed(1)
    n <- 30000
    sampled.new <- sample.dataset(network, n = n)
    set.seed(1) # same seed won't produce identical results (different internals), but that's OK
    sampled.old <- sample.dataset.old(network, n = n)

    rd.new <- raw.data(sampled.new)
    rd.old <- raw.data(sampled.old)

    ns <- node.sizes(network)
    vars <- variables(network)

    for (v in 1:num.nodes(network)) {
        freq.new <- tabulate(rd.new[, v], nbins = ns[v]) / n
        freq.old <- tabulate(rd.old[, v], nbins = ns[v]) / n
        # Both should approximate the true marginal, so they should be close to each other
        expect_equal(freq.new, freq.old,
            tolerance = 0.02,
            info = paste("Marginal mismatch for", vars[v])
        )
    }
})


# --- Performance comparison ---

test_that("vectorized sample.dataset is faster than for-loop version", {
    n <- 5000

    time.new <- system.time({
        for (rep in 1:3) sample.dataset(network, n = n)
    })["elapsed"]

    time.old <- system.time({
        for (rep in 1:3) sample.dataset.old(network, n = n)
    })["elapsed"]

    speedup <- time.old / max(time.new, 0.001)
    message(sprintf(
        "\nPerformance comparison (n=%d, 3 reps):\n  Vectorized: %.3fs\n  For-loop:   %.3fs\n  Speedup:    %.1fx\n",
        n, time.new, time.old, speedup
    ))

    # The vectorized version should be substantially faster
    expect_true(speedup > 5,
        info = sprintf("Expected speedup > 5x, got %.1fx", speedup)
    )
})
