sm_path_subsets <- function(sinks, n.nodes) {
    bitmask <- 2^(seq_len(n.nodes) - 1)
    mask <- 2^n.nodes - 1
    subsets <- integer(n.nodes)

    for (i in seq_len(n.nodes)) {
        subsets[i] <- mask
        sink <- sinks[mask + 1]
        mask <- bitXor(mask, bitmask[sink])
    }

    subsets
}

# ---------------------------------------------------------------------------
# CPDAG (Markov equivalence class) helpers
#
# BDeu is score-equivalent: Markov-equivalent DAGs (same skeleton and same
# v-structures) receive identical scores. Edge orientation inside the optimal
# class is therefore not identifiable from the score, and SM's sink ties are
# exactly this orientation freedom. The well-posed object to certify is the
# CPDAG = (skeleton, v-structures), not the oriented adjacency.
# ---------------------------------------------------------------------------

# undirected skeleton of a DAG adjacency matrix
sm_skeleton <- function(dag) {
    (dag | t(dag)) * 1L
}

# unshielded colliders a -> c <- b with a, b non-adjacent; canonical string keys
sm_vstructures <- function(dag) {
    p <- nrow(dag)
    vs <- character(0)
    for (child in seq_len(p)) {
        pa <- which(dag[, child] == 1L)
        if (length(pa) >= 2) {
            for (a in seq_len(length(pa) - 1)) {
                for (b in seq(a + 1, length(pa))) {
                    x <- pa[a]
                    y <- pa[b]
                    if (dag[x, y] == 0L && dag[y, x] == 0L) {
                        vs <- c(vs, paste0(min(x, y), "-", child, "-", max(x, y)))
                    }
                }
            }
        }
    }
    sort(unique(vs))
}

# CPDAG identity: skeleton + v-structure set
sm_cpdag_key <- function(dag) {
    paste0(paste(sm_skeleton(dag), collapse = ""), "|", paste(sm_vstructures(dag), collapse = ","))
}

# Enumerate score-optimal DAGs by branching on EXACT sink ties.
# Bounded by call.cap / dag.cap; sets truncated = TRUE if either is hit.
#
# IMPORTANT - this branches on sink ties ONLY. At each (node w, subset U) it takes
# the single parent set bps[w, U+1], i.e. whatever fbp's tie-break chose. So when
# two parent sets score EXACTLY equally, the DAGs using the runner-up are never
# produced and the returned list is a strict subset of the optimal class.
#
# That is deliberate and sound, but only because of a coupling with condition 3 of
# sm_cpdag_certificate. Condition 3 requires a STRICT coupled parent margin
# (sm_coupled_parent_witness flags fragile when worst_margin <= 0), which means:
#
#   * an exact parent tie whose coupled term is <= 0 gives margin <= 0, is flagged
#     fragile, and the dataset is not certified - so the un-enumerated branch never
#     matters;
#   * an exact parent tie whose coupled term is > 0 gives margin > 0, and the
#     margin is a valid lower bound on the true post-removal margin, so EVERY
#     removal breaks the tie strictly in P*'s favour. The runner-up is optimal on
#     D but not on D \ {o}, so it cannot appear in any post-removal reconstruction
#     and enumerating P* alone was exactly right.
#
# Either way the DAG the DP rebuilds on D \ {o} is one of the DAGs enumerated here
# (its sinks are D-winners by condition 2, its parent sets are P* by condition 3),
# and condition 1 checked that they all share one CPDAG.
#
# Consequence: do NOT relax condition 3's comparison to a tolerant form such as
# worst_margin < -eps. That would let an exact parent tie (margin exactly 0) pass,
# while this enumeration still ignored the tied alternative - and the certificate
# could then certify a CPDAG that a single removal legitimately changes.
sm_enumerate_optimal_dags <- function(aflml, bps, F, n.nodes, eps = 1e-9,
                                      dag.cap = 20000L, call.cap = 500000L) {
    bit <- 2^(seq_len(n.nodes) - 1)
    full <- 2^n.nodes - 1
    dags <- list()
    seen <- new.env(hash = TRUE)
    # Exact pruning of the recursion. What rec(S, parents) enumerates from here depends
    # only on the remaining set S, and each DAG it emits is the parents fixed so far plus
    # one completion over S. Two branches reaching the same (S, fixed parents) - e.g. the
    # same tied nodes removed in a different order - therefore emit identical DAGs, so
    # the second visit is skipped. Without this the recursion walks every optimal
    # ordering (every linear extension of every optimal DAG), which is combinatorial in
    # the number of isolated or weakly connected nodes and hit call.cap at p = 10 on
    # sparse graphs with only a handful of distinct optima.
    states <- new.env(hash = TRUE)
    truncated <- FALSE
    calls <- 0L

    rec <- function(S, parents) {
        if (truncated) {
            return(invisible())
        }
        skey <- paste(S, paste(vapply(parents, function(x) if (is.null(x)) 0 else x, numeric(1)),
                               collapse = ","), sep = "|")
        if (!is.null(states[[skey]])) {
            return(invisible())
        }
        assign(skey, TRUE, envir = states)
        calls <<- calls + 1L
        if (calls > call.cap) {
            truncated <<- TRUE
            return(invisible())
        }
        if (S == 0) {
            dag <- matrix(0L, n.nodes, n.nodes)
            for (v in seq_len(n.nodes)) {
                ps <- parents[[v]]
                if (!is.null(ps) && ps > 1) {
                    pset <- which(bitAnd(ps - 1, bit) > 0)
                    dag[pset, v] <- 1L
                }
            }
            key <- paste(dag, collapse = "")
            if (is.null(seen[[key]])) {
                assign(key, TRUE, envir = seen)
                dags[[length(dags) + 1]] <<- dag
                if (length(dags) > dag.cap) {
                    truncated <<- TRUE
                }
            }
            return(invisible())
        }
        nodes <- which(bitAnd(S, bit) > 0)
        A <- vapply(nodes, function(i) {
            U <- bitXor(S, bit[i])
            F[U + 1] + aflml[i, bps[i, U + 1]]
        }, numeric(1))
        winners <- nodes[abs(A - max(A)) < eps]
        for (w in winners) {
            if (truncated) {
                return(invisible())
            }
            U <- bitXor(S, bit[w])
            pr <- parents
            pr[[w]] <- bps[w, U + 1]
            rec(U, pr)
        }
    }

    rec(full, vector("list", n.nodes))
    list(dags = dags, truncated = truncated)
}

# Visit every subset reachable as an SM state along ANY optimal-order branch.
# Each subset is expanded once; records the tied winner group at each subset.
sm_optimal_visited <- function(aflml, bps, F, n.nodes, eps = 1e-9) {
    bit <- 2^(seq_len(n.nodes) - 1)
    full <- 2^n.nodes - 1
    visited <- logical(2^n.nodes)
    winfo <- vector("list", 2^n.nodes)
    order_subsets <- integer(0)
    stack <- full
    visited[full + 1] <- TRUE

    while (length(stack) > 0) {
        S <- stack[length(stack)]
        stack <- stack[-length(stack)]
        order_subsets <- c(order_subsets, S)
        nodes <- which(bitAnd(S, bit) > 0)
        A <- vapply(nodes, function(i) {
            U <- bitXor(S, bit[i])
            F[U + 1] + aflml[i, bps[i, U + 1]]
        }, numeric(1))
        winners <- nodes[abs(A - max(A)) < eps]
        winfo[[S + 1]] <- list(nodes = nodes, A = A, winners = winners)
        for (w in winners) {
            U <- bitXor(S, bit[w])
            if (U > 0 && !visited[U + 1]) {
                visited[U + 1] <- TRUE
                stack <- c(stack, U)
            }
        }
    }
    list(subsets = order_subsets, winfo = winfo)
}

# Directional (asymmetric-tight) parent-set fragility test for family (i, U).
# Returns a witness list if the parent choice is not certified, else NULL.
sm_asym_parent_witness <- function(i_star, U_mask, mask, aflml, bps, ifm,
                                   lower_delta, upper_delta) {
    u_idx <- U_mask + 1
    p_star_idx <- bps[i_star, u_idx]
    best_score <- aflml[i_star, p_star_idx]
    if (!is.finite(best_score)) {
        return(NULL)
    }
    max_risk <- -Inf
    p_risk <- NA_integer_
    p_mask <- U_mask
    repeat {
        p_idx <- p_mask + 1
        if (ifm[i_star, p_idx] && p_idx != p_star_idx) {
            risk <- aflml[i_star, p_idx] + upper_delta[i_star, p_idx]
            if (risk > max_risk) {
                max_risk <- risk
                p_risk <- p_idx
            }
        }
        if (p_mask == 0) {
            break
        }
        p_mask <- bitAnd(p_mask - 1, U_mask)
    }
    if (!is.finite(max_risk)) {
        return(NULL)
    }
    lhs <- best_score + lower_delta[i_star, p_star_idx]
    rhs <- max_risk
    if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
        return(list(
            subset = mask,
            node = i_star,
            best_parent_set = p_star_idx,
            competitor_parent_set = p_risk,
            s_best = best_score,
            s_competitor = aflml[i_star, p_risk],
            lower_best = lower_delta[i_star, p_star_idx],
            upper_competitor = upper_delta[i_star, p_risk],
            lhs = lhs,
            rhs = rhs,
            parent_gap = best_score - aflml[i_star, p_risk],
            required_gap = upper_delta[i_star, p_risk] - lower_delta[i_star, p_star_idx],
            type = "parent"
        ))
    }
    NULL
}

# Coupled (same-child-value) parent-set fragility test for family (i, U).
#
# A removed row fixes the child value k across all parent sets of node i, so the
# winner and every competitor must be evaluated at the SAME k. Using per-child
# deltas (lower_k / upper_k) gives, for competitor P,
#   margin(P) = s(P*) - s(P) + min_k [ lower_k(P*,k) - upper_k(P,k) ],
# a valid lower bound on the true leave-one-out margin that is strictly tighter
# than the decoupled min/max bound. The parent choice is certified iff every
# competitor has margin > 0. Returns a witness list if not certified, else NULL.
sm_coupled_parent_witness <- function(i_star, U_mask, mask, aflml, bps, ifm,
                                      node.sizes, lower_k, upper_k,
                                      lower_delta, upper_delta) {
    u_idx <- U_mask + 1
    p_star_idx <- bps[i_star, u_idx]
    best_score <- aflml[i_star, p_star_idx]
    if (!is.finite(best_score)) {
        return(NULL)
    }
    r <- node.sizes[i_star]
    low_star <- lower_k[i_star, p_star_idx, seq_len(r)]

    worst_margin <- Inf
    p_risk <- NA_integer_
    coupled_at_risk <- NA_real_

    p_mask <- U_mask
    repeat {
        p_idx <- p_mask + 1
        if (ifm[i_star, p_idx] && p_idx != p_star_idx) {
            up_P <- upper_k[i_star, p_idx, seq_len(r)]
            valid <- is.finite(low_star) & is.finite(up_P)
            if (any(valid)) {
                coupled <- min(low_star[valid] - up_P[valid])
            } else {
                # no shared child value: fall back to the decoupled bound (safe)
                coupled <- lower_delta[i_star, p_star_idx] - upper_delta[i_star, p_idx]
            }
            # NB the fragile test below is `worst_margin <= 0`, i.e. STRICT
            # certification. sm_enumerate_optimal_dags relies on that: see the
            # note above it before weakening this to a tolerance.
            margin <- best_score - aflml[i_star, p_idx] + coupled
            if (margin < worst_margin) {
                worst_margin <- margin
                p_risk <- p_idx
                coupled_at_risk <- coupled
            }
        }
        if (p_mask == 0) {
            break
        }
        p_mask <- bitAnd(p_mask - 1, U_mask)
    }

    # Always return the decision record (with its margin and a `fragile` flag), so
    # the diagnostic layer can read margins even for certified decisions. The caller
    # records only fragile ones into parent_fragile. NULL means "no decision here"
    # (e.g. the winning family is infeasible or has no competitor).
    if (!is.finite(worst_margin)) {
        return(NULL)
    }
    list(
        subset = mask,
        node = i_star,
        best_parent_set = p_star_idx,
        competitor_parent_set = p_risk,
        s_best = best_score,
        s_competitor = aflml[i_star, p_risk],
        coupled_term = coupled_at_risk,
        margin = worst_margin,
        fragile = worst_margin <= 0,
        type = "parent"
    )
}

# Decoupled (independent min/max) parent margin for family (i, U): the original
# worst-case bound margin(P) = s(P*) - s(P) + lower_delta(P*) - upper_delta(P),
# worst over competitors. Used only by the diagnostic layer to quantify how much
# the coupled-by-k margin tightens over the decoupled one.
sm_decoupled_parent_margin <- function(i_star, U_mask, aflml, bps, ifm,
                                       lower_delta, upper_delta) {
    u_idx <- U_mask + 1
    p_star_idx <- bps[i_star, u_idx]
    s_star <- aflml[i_star, p_star_idx]
    if (!is.finite(s_star)) return(NA_real_)
    worst <- Inf
    p_mask <- U_mask
    repeat {
        p_idx <- p_mask + 1
        if (ifm[i_star, p_idx] && p_idx != p_star_idx) {
            m <- s_star - aflml[i_star, p_idx] +
                lower_delta[i_star, p_star_idx] - upper_delta[i_star, p_idx]
            if (m < worst) worst <- m
        }
        if (p_mask == 0) break
        p_mask <- bitAnd(p_mask - 1, U_mask)
    }
    if (is.finite(worst)) worst else NA_real_
}

# Sound CPDAG-stability certificate (see stability_findings_and_redesign.md, 3.1).
# Uses asymmetric-tight bounds Lg/Ug/LF/UF plus the coupled local delta bounds.
sm_cpdag_certificate <- function(aflml, bps, F, ifm, lower_delta, upper_delta,
                                 lower_k, upper_k, node.sizes,
                                 Lg, Ug, LF, UF, dag, n.nodes, eps = 1e-9,
                                 diagnostics = FALSE) {
    bit <- 2^(seq_len(n.nodes) - 1)
    base_key <- sm_cpdag_key(dag)

    # Condition 1: the optimal Markov equivalence class must be a single CPDAG.
    enum <- sm_enumerate_optimal_dags(aflml, bps, F, n.nodes, eps)
    cpdag_keys <- vapply(enum$dags, sm_cpdag_key, character(1))
    optimal_cpdag_unique <- (!enum$truncated) && length(cpdag_keys) > 0 &&
        all(cpdag_keys == base_key)
    skel_keys <- unique(vapply(enum$dags, function(d) paste(sm_skeleton(d), collapse = ""), character(1)))
    vs_keys <- unique(vapply(enum$dags, function(d) paste(sm_vstructures(d), collapse = ","), character(1)))

    # Conditions 2 & 3 over the union of all optimal-order subsets.
    vis <- sm_optimal_visited(aflml, bps, F, n.nodes, eps)
    sink_fragile <- list()
    parent_fragile <- list()
    # Diagnostic records (only populated when diagnostics = TRUE). These reuse the
    # margins the certificate already computes, so the extra cost is O(1) per
    # decision; the default path stores nothing.
    sink_diag <- if (diagnostics) list() else NULL
    parent_diag <- if (diagnostics) list() else NULL

    for (S in vis$subsets) {
        info <- vis$winfo[[S + 1]]
        nodes <- info$nodes
        winners <- info$winners

        # Condition 2: no strict competitor (original gap > 0) may be promoted
        # above the whole tied winner group after any single removal.
        if (length(nodes) >= 2 && length(winners) < length(nodes)) {
            lhs <- -Inf
            A_wbest <- NA_real_
            LA_wbest <- NA_real_
            for (w in winners) {
                U <- bitXor(S, bit[w])
                A_w <- F[U + 1] + aflml[w, bps[w, U + 1]]
                LA_w <- LF[U + 1] + Lg[w, U + 1]
                val <- A_w + LA_w
                if (val > lhs) { lhs <- val; A_wbest <- A_w; LA_wbest <- LA_w }
            }
            strict <- setdiff(nodes, winners)
            rhs <- -Inf
            j_risk <- NA_integer_
            A_jrisk <- NA_real_; UF_r <- NA_real_; Ug_r <- NA_real_
            for (j in strict) {
                U <- bitXor(S, bit[j])
                A_j <- F[U + 1] + aflml[j, bps[j, U + 1]]
                UF_j <- UF[U + 1]; Ug_j <- Ug[j, U + 1]
                val <- A_j + UF_j + Ug_j
                if (val > rhs) {
                    rhs <- val; j_risk <- j
                    A_jrisk <- A_j; UF_r <- UF_j; Ug_r <- Ug_j
                }
            }
            if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
                sink_fragile[[length(sink_fragile) + 1]] <- list(
                    subset = S,
                    winners = winners,
                    competitor = j_risk,
                    lhs = lhs,
                    rhs = rhs,
                    gap = lhs - rhs,
                    type = "sink_cpdag"
                )
            }
            # sink-margin anatomy: margin = base_gap + LA_win - (UF + Ug)
            if (diagnostics) {
                sink_diag[[length(sink_diag) + 1]] <- data.frame(
                    subset = S, competitor = j_risk,
                    margin = lhs - rhs,
                    base_gap = A_wbest - A_jrisk,
                    UF = UF_r, Ug = Ug_r, LA_win = LA_wbest
                )
            }
        }

        # Condition 3: family/parent decision stable for each optimal winner
        # (coupled same-child-value margin).
        for (w in winners) {
            U <- bitXor(S, bit[w])
            wt <- sm_coupled_parent_witness(
                w, U, S, aflml, bps, ifm, node.sizes,
                lower_k, upper_k, lower_delta, upper_delta
            )
            if (!is.null(wt)) {
                if (isTRUE(wt$fragile)) {
                    parent_fragile[[length(parent_fragile) + 1]] <- wt
                }
                if (diagnostics) {
                    dec <- sm_decoupled_parent_margin(
                        w, U, aflml, bps, ifm, lower_delta, upper_delta
                    )
                    pstar_size <- sum(bitAnd(bps[w, U + 1] - 1, bit) > 0)
                    parent_diag[[length(parent_diag) + 1]] <- data.frame(
                        subset = S, node = w, pstar_size = pstar_size,
                        margin_coupled = wt$margin, margin_decoupled = dec
                    )
                }
            }
        }
    }

    stable_cpdag <- optimal_cpdag_unique &&
        length(sink_fragile) == 0 && length(parent_fragile) == 0

    out <- list(
        stable_cpdag = stable_cpdag,
        optimal_cpdag_unique = optimal_cpdag_unique,
        skeleton_unique_among_optima = length(skel_keys) == 1,
        vstructures_unique_among_optima = length(vs_keys) == 1,
        n_optimal_dags = length(enum$dags),
        enumeration_truncated = enum$truncated,
        sink_fragile = sink_fragile,
        parent_fragile = parent_fragile,
        witnesses = c(sink_fragile, parent_fragile),
        base_cpdag_key = base_key
    )

    if (diagnostics) {
        sink_df <- if (length(sink_diag)) do.call(rbind, sink_diag) else NULL
        parent_df <- if (length(parent_diag)) do.call(rbind, parent_diag) else NULL
        min_sink <- if (!is.null(sink_df)) min(sink_df$margin) else Inf
        min_parent <- if (!is.null(parent_df)) min(parent_df$margin_coupled) else Inf
        binding <- if (min_sink == Inf && min_parent == Inf) "none" else
            if (min_sink <= min_parent) "sink" else "parent"
        out$diagnostics <- list(
            sink = sink_df,
            parent = parent_df,
            min_sink_margin = min_sink,
            min_parent_margin = min_parent,
            binding = binding
        )
    }
    out
}

# Recover the exact SM-optimal DAG from a family-score table (BDeu, encoded by
# the same importance-family mask). Mirrors the core of sm().
sm_optimal_dag_from_aflml <- function(aflml, ifm) {
    bps <- .Call("bnstruct_fbp", aflml = aflml, PACKAGE = "bnstruct")
    sinks <- .Call("bnstruct_fbs", bps = bps, aflml = aflml, PACKAGE = "bnstruct")
    order <- find.best.ordering(sinks)
    parents <- find.best.network(order, bps)
    parents.to.dag(parents)
}

# Exact incremental leave-one-out CPDAG certificate. For every observation o, the
# optimal DAG on D \ {o} is recovered by the C routine bnstruct_loo_parents, which
# shifts the full-data family scores by the closed-form delta of the single cell
# row o occupies and re-runs only the SM dynamic program (no re-counting). Its
# CPDAG (skeleton + v-structures) is then compared to the base CPDAG. Necessary
# and sufficient (no bound slack); also returns the structurally influential
# observations and what each of them changes.
sm_loo_exact_cpdag <- function(data, node.sizes, ess, ifm, dag, aflml,
                               scoring.func = 0L) {
    n <- nrow(data)
    p <- ncol(data)
    base_key <- sm_cpdag_key(dag)
    base_skel <- sm_skeleton(dag)
    base_vs <- sm_vstructures(dag)

    # parents_mat[o, i] = 0-indexed parent-set bitmask of node i in the optimum for D \ {o}
    parents_mat <- .Call(
        "bnstruct_loo_parents_score",
        data, node.sizes, ifm, ess, aflml, as.integer(scoring.func),
        PACKAGE = "bnstruct"
    )
    bit <- 2^(seq_len(p) - 1)

    # Base parent-set bitmask per node. If a row's parent vector equals this, the
    # DAG is byte-identical, hence the CPDAG is identical: skip the (comparatively
    # expensive) v-structure / CPDAG-key computation. Only rows whose parent vector
    # actually changed are examined further. This is exact, not a heuristic prune.
    base_pm <- vapply(seq_len(p), function(i) {
        pa <- which(dag[, i] == 1L)
        if (length(pa) == 0) 0L else as.integer(sum(bit[pa]))
    }, integer(1))
    changed_rows <- which(!apply(parents_mat, 1, function(r) all(r == base_pm)))

    influential <- list()
    n_skeleton_changes <- 0L
    n_vstructure_changes <- 0L

    for (o in changed_rows) {
        dag_o <- matrix(0L, p, p)
        for (i in seq_len(p)) {
            pm <- parents_mat[o, i]
            if (pm > 0) {
                dag_o[which(bitAnd(pm, bit) > 0), i] <- 1L
            }
        }

        if (sm_cpdag_key(dag_o) != base_key) {
            skel_o <- sm_skeleton(dag_o)
            vs_o <- sm_vstructures(dag_o)
            skel_changed <- !identical(skel_o, base_skel)
            vs_changed <- !identical(vs_o, base_vs)
            if (skel_changed) n_skeleton_changes <- n_skeleton_changes + 1L
            if (vs_changed) n_vstructure_changes <- n_vstructure_changes + 1L
            influential[[length(influential) + 1]] <- list(
                obs = o,
                skeleton_changed = skel_changed,
                vstructures_changed = vs_changed,
                added_vstructures = setdiff(vs_o, base_vs),
                removed_vstructures = setdiff(base_vs, vs_o)
            )
        }
    }

    list(
        stable_cpdag = length(influential) == 0,
        target = "cpdag",
        mode = "exact",
        n_obs = n,
        n_influential = length(influential),
        n_skeleton_changes = n_skeleton_changes,
        n_vstructure_changes = n_vstructure_changes,
        influential = influential,
        base_cpdag_key = base_key
    )
}

sm_stability_bdeu <- function(
    data,
    node.sizes,
    ess,
    aflml,
    bps,
    sinks,
    F,
    ifm,
    dag = NULL,
    target = c("cpdag", "dag"),
    mode = c("bound", "exact"),
    bound.type = c("cons", "tight", "asym_cons", "asym_tight"),
    scope = c("path", "all"),
    diagnostics = FALSE,
    scoring.func = 0L,
    bound.supplied = FALSE,
    scope.supplied = FALSE) {
    # scoring.func: 0 = BDeu (ess used), 1 = AIC, 2 = BIC (ess ignored). Every
    # certificate below only needs (a) the per-family removal-delta bounds and
    # (b) delta >= 0 and score equivalence, all of which hold for the three
    # scores; aflml must have been computed with the same scoring.func.
    #
    # bound.supplied / scope.supplied say whether the caller set bound.type /
    # scope explicitly (sm() derives them with missing()). Both arguments apply
    # to target = "dag" only; accepting them silently under target = "cpdag"
    # would suggest a control the CPDAG certificate does not have.
    target <- match.arg(target)
    mode <- match.arg(mode)
    bound.type <- match.arg(bound.type)
    scope <- match.arg(scope)
    scoring.func <- as.integer(scoring.func)

    if (target == "cpdag") {
        if (bound.supplied) {
            stop(
                "stability.bound applies only to stability.target = 'dag'.\n",
                "  The CPDAG certificate always uses the asymmetric-tight bounds ",
                "(Lg/Ug via bnstruct_compute_lg_ug_tight, LF/UF via bnstruct_compute_lf_uf) ",
                "together with the coupled per-child-value parent margin; the symmetric ",
                "'cons'/'tight' variants have no CPDAG counterpart.\n",
                "  Drop stability.bound, or set stability.target = 'dag'."
            )
        }
        if (scope.supplied) {
            stop(
                "stability.scope applies only to stability.target = 'dag'.\n",
                "  The CPDAG certificate derives its own set of subsets - the union over ",
                "every optimal-order branch (sm_optimal_visited) - which is neither the ",
                "single backtracking path ('path') nor all 2^p subsets ('all').\n",
                "  Drop stability.scope, or set stability.target = 'dag'."
            )
        }
    } else {
        # scope = "path" walks only the subsets on the backtracking path
        # S_0 = V, S_{t+1} = S_t \ {Last(S_t)}. That is sound for the oriented-DAG
        # target: certifying the sink at every path subset means the perturbed DP
        # follows the same path, so no off-path subset is ever visited. "all"
        # would check every subset instead - a strictly stronger statement about
        # the whole DP table, not needed for DAG stability.
        if (scope != "path") {
            stop("stability.scope = 'all' not yet implemented.")
        }

        # Exact leave-one-out is implemented for the CPDAG target only. Silently
        # falling back to the bound certificate would downgrade a necessary-and-
        # sufficient answer to a merely sufficient one, so refuse instead.
        if (mode == "exact") {
            stop("stability.mode = 'exact' is implemented only for stability.target = 'cpdag'.")
        }
    }

    n.nodes <- ncol(data)

    bounds <- .Call(
        "bnstruct_score_loo_bounds",
        data, node.sizes, ifm, ess, scoring.func,
        PACKAGE = "bnstruct"
    )

    lower_delta <- bounds$lower
    upper_delta <- bounds$upper

    # CPDAG (Markov equivalence class) certificate. Certifies invariance of the
    # skeleton and v-structures rather than the (unidentifiable) oriented DAG.
    if (target == "cpdag") {
        if (is.null(dag)) {
            stop("target = 'cpdag' requires the learned dag.")
        }

        # Exact incremental leave-one-out: re-run the SM optimum on D \ {o} for
        # every row and compare CPDAGs. Necessary and sufficient (no bound slack)
        # and yields the list of structurally influential observations.
        if (mode == "exact") {
            cert <- sm_loo_exact_cpdag(data, node.sizes, ess, ifm, dag, aflml,
                                       scoring.func = scoring.func)
            # exact mode examines every observation; no subset restriction applies
            cert$scope <- "all_observations"
            return(cert)
        }

        g_bounds <- .Call(
            "bnstruct_compute_lg_ug_tight",
            aflml, bps, lower_delta, upper_delta, ifm,
            PACKAGE = "bnstruct"
        )
        f_bounds <- .Call(
            "bnstruct_compute_lf_uf",
            F, sinks, bps, aflml, g_bounds$lower, g_bounds$upper,
            PACKAGE = "bnstruct"
        )
        cert <- sm_cpdag_certificate(
            aflml, bps, F, ifm, lower_delta, upper_delta,
            bounds$lower_k, bounds$upper_k, node.sizes,
            g_bounds$lower, g_bounds$upper, f_bounds$lower, f_bounds$upper,
            dag, n.nodes, diagnostics = diagnostics
        )
        cert$target <- "cpdag"
        cert$bound_type <- "asym_tight"
        # the certificate checks every subset reachable on any optimal-order
        # branch, not the single backtracking path
        cert$scope <- "optimal_visited"
        return(cert)
    }

    is_asym <- bound.type %in% c("asym_cons", "asym_tight")
    is_tight <- bound.type %in% c("tight", "asym_tight")

    if (!is_asym) {
        if (bound.type == "cons") {
            Bg <- .Call(
                "bnstruct_compute_bg_cons",
                lower_delta, upper_delta, ifm,
                PACKAGE = "bnstruct"
            )
        } else {
            Bg <- .Call(
                "bnstruct_compute_bg_tight",
                aflml, bps, lower_delta, upper_delta, ifm,
                PACKAGE = "bnstruct"
            )
        }

        BF <- .Call(
            "bnstruct_compute_bf",
            Bg,
            PACKAGE = "bnstruct"
        )

        blocal <- pmax(abs(lower_delta), abs(upper_delta))
    } else {
        if (bound.type == "asym_cons") {
            g_bounds <- .Call(
                "bnstruct_compute_lg_ug_cons",
                bps, lower_delta, upper_delta, ifm,
                PACKAGE = "bnstruct"
            )
        } else {
            g_bounds <- .Call(
                "bnstruct_compute_lg_ug_tight",
                aflml, bps, lower_delta, upper_delta, ifm,
                PACKAGE = "bnstruct"
            )
        }

        Lg <- g_bounds$lower
        Ug <- g_bounds$upper

        f_bounds <- .Call(
            "bnstruct_compute_lf_uf",
            F, sinks, bps, aflml, Lg, Ug,
            PACKAGE = "bnstruct"
        )

        LF <- f_bounds$lower
        UF <- f_bounds$upper
    }

    n.nodes <- ncol(data)
    bitmask <- 2^(seq_len(n.nodes) - 1)

    subsets <- sm_path_subsets(sinks, n.nodes)

    sink_failures <- list()
    parent_failures <- list()

    for (mask in subsets) {
        nodes <- which(bitAnd(mask, bitmask) > 0)
        if (length(nodes) == 0) {
            next
        }

        i_star <- sinks[mask + 1]

        if (length(nodes) >= 2) {
            A_vals <- numeric(length(nodes))

            if (is_asym) {
                LA_vals <- numeric(length(nodes))
                UA_vals <- numeric(length(nodes))
            } else {
                BA_vals <- numeric(length(nodes))
            }

            for (idx in seq_along(nodes)) {
                i <- nodes[idx]
                U_mask <- bitXor(mask, bitmask[i])
                u_idx <- U_mask + 1
                A_vals[idx] <- F[u_idx] + aflml[i, bps[i, u_idx]]

                if (is_asym) {
                    LA_vals[idx] <- LF[u_idx] + Lg[i, u_idx]
                    UA_vals[idx] <- UF[u_idx] + Ug[i, u_idx]
                } else {
                    BA_vals[idx] <- BF[u_idx] + Bg[i, u_idx]
                }
            }

            idx_star <- which(nodes == i_star)
            if (length(idx_star) == 1) {
                comp_idx <- setdiff(seq_along(nodes), idx_star)
                if (length(comp_idx) > 0) {
                    if (is_asym) {
                        risk_vals <- A_vals[comp_idx] + UA_vals[comp_idx]
                        if (any(is.finite(risk_vals))) {
                            risk_idx <- comp_idx[which.max(risk_vals)]

                            A_star <- A_vals[idx_star]
                            LA_star <- LA_vals[idx_star]
                            UA_star <- UA_vals[idx_star]
                            lhs <- A_star + LA_star
                            rhs <- A_vals[risk_idx] + UA_vals[risk_idx]

                            if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
                                sink_failures[[length(sink_failures) + 1]] <- list(
                                    subset = mask,
                                    chosen = i_star,
                                    competitor = nodes[risk_idx],
                                    A_chosen = A_star,
                                    A_competitor = A_vals[risk_idx],
                                    LA_chosen = LA_star,
                                    UA_chosen = UA_star,
                                    LA_competitor = LA_vals[risk_idx],
                                    UA_competitor = UA_vals[risk_idx],
                                    lhs = lhs,
                                    rhs = rhs,
                                    gap = A_star - A_vals[risk_idx],
                                    required_gap = UA_vals[risk_idx] - LA_star,
                                    type = "sink"
                                )
                            }
                        }
                    } else if (is_tight) {
                        risk_vals <- A_vals[comp_idx] + BA_vals[comp_idx]
                        risk_idx <- comp_idx[which.max(risk_vals)]

                        A_star <- A_vals[idx_star]
                        BA_star <- BA_vals[idx_star]
                        lhs <- A_star - BA_star
                        rhs <- A_vals[risk_idx] + BA_vals[risk_idx]

                        if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
                            sink_failures[[length(sink_failures) + 1]] <- list(
                                subset = mask,
                                chosen = i_star,
                                competitor = nodes[risk_idx],
                                A_chosen = A_star,
                                A_competitor = A_vals[risk_idx],
                                BA_chosen = BA_star,
                                BA_competitor = BA_vals[risk_idx],
                                lhs = lhs,
                                rhs = rhs,
                                gap = A_star - A_vals[risk_idx],
                                bound = BA_star + BA_vals[risk_idx],
                                type = "sink"
                            )
                        }
                    } else {
                        comp_scores <- A_vals[comp_idx]
                        if (any(is.finite(comp_scores))) {
                            best_comp_idx <- comp_idx[which.max(comp_scores)]
                            A_star <- A_vals[idx_star]
                            A_second <- A_vals[best_comp_idx]
                            BA_star <- BA_vals[idx_star]
                            BA_second <- BA_vals[best_comp_idx]
                            gap <- A_star - A_second
                            bound <- BA_star + BA_second

                            if (!is.finite(gap) || !is.finite(bound) || gap <= bound) {
                                sink_failures[[length(sink_failures) + 1]] <- list(
                                    subset = mask,
                                    chosen = i_star,
                                    runner_up = nodes[best_comp_idx],
                                    gap = gap,
                                    bound = bound,
                                    type = "sink"
                                )
                            }
                        }
                    }
                }
            }
        }

        U_mask <- bitXor(mask, bitmask[i_star])
        u_idx <- U_mask + 1
        p_star_idx <- bps[i_star, u_idx]
        best_score <- aflml[i_star, p_star_idx]

        if (is.finite(best_score)) {
            if (is_asym) {
                max_risk <- -Inf
                p_risk <- NA_integer_

                p_mask <- U_mask
                repeat {
                    p_idx <- p_mask + 1
                    if (ifm[i_star, p_idx] && p_idx != p_star_idx) {
                        risk <- aflml[i_star, p_idx] + upper_delta[i_star, p_idx]
                        if (risk > max_risk) {
                            max_risk <- risk
                            p_risk <- p_idx
                        }
                    }

                    if (p_mask == 0) {
                        break
                    }
                    p_mask <- bitAnd(p_mask - 1, U_mask)
                }

                if (is.finite(max_risk)) {
                    lhs <- best_score + lower_delta[i_star, p_star_idx]
                    rhs <- max_risk
                    parent_gap <- best_score - aflml[i_star, p_risk]
                    required_gap <- upper_delta[i_star, p_risk] - lower_delta[i_star, p_star_idx]

                    if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
                        parent_failures[[length(parent_failures) + 1]] <- list(
                            subset = mask,
                            node = i_star,
                            best_parent_set = p_star_idx,
                            competitor_parent_set = p_risk,
                            s_best = best_score,
                            s_competitor = aflml[i_star, p_risk],
                            lower_best = lower_delta[i_star, p_star_idx],
                            upper_competitor = upper_delta[i_star, p_risk],
                            lhs = lhs,
                            rhs = rhs,
                            parent_gap = parent_gap,
                            required_gap = required_gap,
                            type = "parent"
                        )
                    }
                }
            } else if (is_tight) {
                max_risk <- -Inf
                p_risk <- NA_integer_

                p_mask <- U_mask
                repeat {
                    p_idx <- p_mask + 1
                    if (ifm[i_star, p_idx] && p_idx != p_star_idx) {
                        risk <- aflml[i_star, p_idx] + upper_delta[i_star, p_idx]
                        if (risk > max_risk) {
                            max_risk <- risk
                            p_risk <- p_idx
                        }
                    }

                    if (p_mask == 0) {
                        break
                    }
                    p_mask <- bitAnd(p_mask - 1, U_mask)
                }

                if (is.finite(max_risk)) {
                    lhs <- best_score + lower_delta[i_star, p_star_idx]
                    rhs <- max_risk
                    parent_gap <- best_score - aflml[i_star, p_risk]
                    parent_bound <- upper_delta[i_star, p_risk] - lower_delta[i_star, p_star_idx]

                    if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
                        parent_failures[[length(parent_failures) + 1]] <- list(
                            subset = mask,
                            node = i_star,
                            best_parent_set = p_star_idx,
                            competitor_parent_set = p_risk,
                            lhs = lhs,
                            rhs = rhs,
                            parent_gap = parent_gap,
                            bound = parent_bound,
                            type = "parent"
                        )
                    }
                }
            } else {
                max_other_score <- -Inf
                max_other_blocal <- -Inf
                blocal_star <- blocal[i_star, p_star_idx]

                p_mask <- U_mask
                repeat {
                    p_idx <- p_mask + 1
                    if (ifm[i_star, p_idx]) {
                        score <- aflml[i_star, p_idx]
                        bl <- blocal[i_star, p_idx]
                        if (p_idx != p_star_idx) {
                            if (score > max_other_score) {
                                max_other_score <- score
                            }
                            if (bl > max_other_blocal) {
                                max_other_blocal <- bl
                            }
                        }
                    }

                    if (p_mask == 0) {
                        break
                    }
                    p_mask <- bitAnd(p_mask - 1, U_mask)
                }

                if (is.finite(max_other_score)) {
                    parent_gap <- best_score - max_other_score
                    parent_bound <- blocal_star + max_other_blocal

                    if (!is.finite(parent_gap) || !is.finite(parent_bound) || parent_gap <= parent_bound) {
                        parent_failures[[length(parent_failures) + 1]] <- list(
                            subset = mask,
                            node = i_star,
                            best_parent_set = p_star_idx,
                            parent_gap = parent_gap,
                            bound = parent_bound,
                            type = "parent"
                        )
                    }
                }
            }
        }
    }

    stable_order <- length(sink_failures) == 0
    stable_dag <- stable_order && length(parent_failures) == 0

    list(
        stable_order = stable_order,
        stable_dag = stable_dag,
        sink_failures = sink_failures,
        parent_failures = parent_failures,
        witnesses = c(sink_failures, parent_failures),
        bound_type = bound.type,
        scope = scope
    )
}
