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

# Enumerate all score-optimal DAGs by branching on EXACT sink ties.
# Bounded by call.cap / dag.cap; sets truncated = TRUE if either is hit.
sm_enumerate_optimal_dags <- function(aflml, bps, F, n.nodes, eps = 1e-9,
                                      dag.cap = 20000L, call.cap = 500000L) {
    bit <- 2^(seq_len(n.nodes) - 1)
    full <- 2^n.nodes - 1
    dags <- list()
    seen <- new.env(hash = TRUE)
    truncated <- FALSE
    calls <- 0L

    rec <- function(S, parents) {
        if (truncated) {
            return(invisible())
        }
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

    if (is.finite(worst_margin) && worst_margin <= 0) {
        return(list(
            subset = mask,
            node = i_star,
            best_parent_set = p_star_idx,
            competitor_parent_set = p_risk,
            s_best = best_score,
            s_competitor = aflml[i_star, p_risk],
            coupled_term = coupled_at_risk,
            margin = worst_margin,
            type = "parent"
        ))
    }
    NULL
}

# Sound CPDAG-stability certificate (see stability_findings_and_redesign.md, 3.1).
# Uses asymmetric-tight bounds Lg/Ug/LF/UF plus the coupled local delta bounds.
sm_cpdag_certificate <- function(aflml, bps, F, ifm, lower_delta, upper_delta,
                                 lower_k, upper_k, node.sizes,
                                 Lg, Ug, LF, UF, dag, n.nodes, eps = 1e-9) {
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

    for (S in vis$subsets) {
        info <- vis$winfo[[S + 1]]
        nodes <- info$nodes
        winners <- info$winners

        # Condition 2: no strict competitor (original gap > 0) may be promoted
        # above the whole tied winner group after any single removal.
        if (length(nodes) >= 2 && length(winners) < length(nodes)) {
            lhs <- -Inf
            for (w in winners) {
                U <- bitXor(S, bit[w])
                A_w <- F[U + 1] + aflml[w, bps[w, U + 1]]
                LA_w <- LF[U + 1] + Lg[w, U + 1]
                val <- A_w + LA_w
                if (val > lhs) lhs <- val
            }
            strict <- setdiff(nodes, winners)
            rhs <- -Inf
            j_risk <- NA_integer_
            for (j in strict) {
                U <- bitXor(S, bit[j])
                A_j <- F[U + 1] + aflml[j, bps[j, U + 1]]
                UA_j <- UF[U + 1] + Ug[j, U + 1]
                val <- A_j + UA_j
                if (val > rhs) {
                    rhs <- val
                    j_risk <- j
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
                parent_fragile[[length(parent_fragile) + 1]] <- wt
            }
        }
    }

    stable_cpdag <- optimal_cpdag_unique &&
        length(sink_fragile) == 0 && length(parent_fragile) == 0

    list(
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
sm_loo_exact_cpdag <- function(data, node.sizes, ess, ifm, dag, aflml) {
    n <- nrow(data)
    p <- ncol(data)
    base_key <- sm_cpdag_key(dag)
    base_skel <- sm_skeleton(dag)
    base_vs <- sm_vstructures(dag)

    # parents_mat[o, i] = 0-indexed parent-set bitmask of node i in the optimum for D \ {o}
    parents_mat <- .Call(
        "bnstruct_loo_parents",
        data, node.sizes, ifm, ess, aflml,
        PACKAGE = "bnstruct"
    )
    bit <- 2^(seq_len(p) - 1)

    influential <- list()
    n_skeleton_changes <- 0L
    n_vstructure_changes <- 0L

    for (o in seq_len(n)) {
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
    scope = c("path", "all")) {
    target <- match.arg(target)
    mode <- match.arg(mode)
    bound.type <- match.arg(bound.type)
    scope <- match.arg(scope)

    if (target == "dag" && scope != "path") {
        stop("stability.scope = 'all' not yet implemented.")
    }

    n.nodes <- ncol(data)

    bounds <- .Call(
        "bnstruct_bdeu_loo_bounds",
        data, node.sizes, ifm, ess,
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
            cert <- sm_loo_exact_cpdag(data, node.sizes, ess, ifm, dag, aflml)
            cert$scope <- scope
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
            dag, n.nodes
        )
        cert$target <- "cpdag"
        cert$bound_type <- "asym_tight"
        cert$scope <- scope
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
